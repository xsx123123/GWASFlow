#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GWASFlow 分阶段基因型质控与群体结构分析（参数手册 §3 契约生产化）。

--mode qc        四阶段独立执行（各自独立 PLINK2 调用与产物，禁止单混合命令）：
                 stage1 初步位点缺失(pre_geno) → stage2 样本缺失(mind)
                 → stage3 最终位点缺失/MAF/MAC → final_stats 最终队列重新统计
--mode structure LD 剪枝 + 剪枝集 PCA（关联分析仍用完整 cohort 标记集）

固定非人类语义：--chr-set {chr_count} no-x no-y no-xy no-mt --nonfounders；
HWE 默认关闭且审计记录；高亲缘材料保留（不自动删除）。

错误码：0=成功；2=格式/配置/输入错误；外部软件非零退出透传。
"""

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

EXIT_OK = 0
EXIT_INPUT_ERROR = 2


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="genotype_qc.py",
        description="GWASFlow 分阶段基因型 QC（--mode qc）与 LD 剪枝/PCA（--mode structure）。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--mode", default="qc", choices=["qc", "structure"],
                        help="执行模式：qc=四阶段质控；structure=剪枝+PCA")
    parser.add_argument("--bfile", required=True, help="输入 BED 前缀（不含扩展名）")
    parser.add_argument("--gate", default=None, help="01.gate/input_gate.json（门控产物）")
    parser.add_argument("--chr-count", type=int, default=9,
                        help="二倍体常染色体对数，1..95")
    parser.add_argument("--pre-geno", type=float, default=0.20, help="初步宽松位点缺失阈值（比例）")
    parser.add_argument("--mind", type=float, default=0.10, help="样本缺失阈值（比例）")
    parser.add_argument("--geno", type=float, default=0.10, help="最终位点缺失阈值（比例）")
    parser.add_argument("--maf", type=float, default=0.05, help="队列级最小 MAF（比例，<=0.5）")
    parser.add_argument("--mac", type=int, default=5, help="队列级最小 MAC（正整数计数）")
    parser.add_argument("--min-samples-after-qc", type=int, default=20,
                        help="QC 后最少样本数（阻断阈值）")
    parser.add_argument("--threads", type=int, default=4, help="PLINK 线程上限")
    parser.add_argument("--memory-mb", type=int, default=8000, help="PLINK --memory 工作内存")
    # structure 模式参数
    parser.add_argument("--prune-window-kb", type=int, default=50, help="剪枝物理窗口（kb）")
    parser.add_argument("--prune-step", type=int, default=5, help="剪枝步长（标记数）")
    parser.add_argument("--prune-r2", type=float, default=0.20, help="剪枝 r2 阈值 (0,1]")
    parser.add_argument("--pca-components", type=int, default=10, help="输出 PC 数")
    parser.add_argument("--out-dir", required=True, help="输出目录（必须为新目录或空目录）")
    return parser.parse_args(argv)


def fail(message, code=EXIT_INPUT_ERROR):
    print(f"[genotype_qc ERROR] {message}", file=sys.stderr)
    sys.exit(code)


def ensure_empty_dir(out_dir):
    path = Path(out_dir).resolve()
    if path.exists():
        if not path.is_dir() or any(path.iterdir()):
            fail(f"输出目录非空或非目录（拒绝覆盖既有产物）: {path}")
    else:
        path.mkdir(parents=True, exist_ok=True)
    return path


def run_plink2(argv, log_path, threads, memory_mb):
    """执行 plink2 命令数组（固定非人类语义参数由调用方组装）；退出码透传。"""
    env = os.environ.copy()
    full = ["plink2", *argv, "--threads", str(threads), "--memory", str(memory_mb)]
    with open(log_path, "a", encoding="utf-8") as handle:
        handle.write(f"\n=== RUN: {' '.join(full)} ===\n")
        handle.flush()
        result = subprocess.run(full, stdout=handle, stderr=subprocess.STDOUT,
                                env=env, check=False)
    return result.returncode


def chr_set_args(chr_count):
    """固定非人类染色体集语义：no-x no-y no-xy no-mt。"""
    return ["--chr-set", str(chr_count), "no-x", "no-y", "no-xy", "no-mt", "--nonfounders"]


def read_bim_ids(bim_path):
    ids = []
    with open(bim_path, encoding="utf-8") as handle:
        for line in handle:
            parts = line.split()
            if len(parts) >= 2:
                ids.append(parts[1])
    return ids


def read_fam_ids(fam_path):
    ids = []
    with open(fam_path, encoding="utf-8") as handle:
        for line in handle:
            parts = line.split()
            if len(parts) >= 2:
                ids.append((parts[0], parts[1]))
    return ids


def parse_plink_log_reasons(log_path):
    """解析 plink2 log 中 *removed due to* 原因行（补充确定性集合对比）。"""
    reasons = []
    if not Path(log_path).exists():
        return reasons
    with open(log_path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if "removed due to" in line.lower() or "removed because" in line.lower():
                reasons.append(line.strip())
    return reasons


def write_tsv(path, header, rows):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)


def stage_stats_row(stage, n_samples_before, n_variants_before,
                    n_samples_after, n_variants_after, params, reasons):
    return {
        "stage": stage,
        "n_samples_before": n_samples_before,
        "n_variants_before": n_variants_before,
        "n_samples_removed": n_samples_before - n_samples_after,
        "n_variants_removed": n_variants_before - n_variants_after,
        "n_samples_after": n_samples_after,
        "n_variants_after": n_variants_after,
        "params": params,
        "reasons": reasons,
    }


STATS_HEADER = [
    "stage", "n_samples_before", "n_variants_before", "n_samples_removed",
    "n_variants_removed", "n_samples_after", "n_variants_after", "params",
]


def run_qc(args, out_dir):
    log_path = out_dir / "plink2.log"
    raw = args.bfile

    n_samples = len(read_fam_ids(f"{raw}.fam"))
    n_variants = len(read_bim_ids(f"{raw}.bim"))

    # ---- stage1: 初步宽松位点缺失 ----
    stage1_prefix = out_dir / "stage1"
    code = run_plink2(
        ["--bfile", raw, "--geno", str(args.pre_geno), *chr_set_args(args.chr_count),
         "--make-bed", "--out", str(stage1_prefix)],
        log_path, args.threads, args.memory_mb,
    )
    if code != 0:
        fail(f"stage1（pre-geno={args.pre_geno}）plink2 失败（退出码 {code}），日志: {log_path}")
    bim_before, bim_after = read_bim_ids(f"{raw}.bim"), read_bim_ids(f"{stage1_prefix}.bim")
    removed_v1 = sorted(set(bim_before) - set(bim_after))
    write_tsv(out_dir / "removed_variants_stage1.tsv",
              ["variant_id"], [[v] for v in removed_v1])
    s1 = stage_stats_row("stage1_pre_geno", n_samples, n_variants,
                         len(read_fam_ids(f"{stage1_prefix}.fam")), len(bim_after),
                         {"geno": args.pre_geno},
                         parse_plink_log_reasons(f"{stage1_prefix}.log"))

    # ---- stage2: 样本缺失 ----
    stage2_prefix = out_dir / "stage2"
    code = run_plink2(
        ["--bfile", str(stage1_prefix), "--mind", str(args.mind), *chr_set_args(args.chr_count),
         "--make-bed", "--out", str(stage2_prefix)],
        log_path, args.threads, args.memory_mb,
    )
    if code != 0:
        fail(f"stage2（mind={args.mind}）plink2 失败（退出码 {code}），日志: {log_path}")
    fam_before, fam_after = read_fam_ids(f"{stage1_prefix}.fam"), read_fam_ids(f"{stage2_prefix}.fam")
    removed_s2 = sorted(set(fam_before) - set(fam_after))
    write_tsv(out_dir / "removed_samples_stage2.tsv",
              ["fid", "iid"], [list(p) for p in removed_s2])
    s2 = stage_stats_row("stage2_mind", s1["n_samples_after"], s1["n_variants_after"],
                         len(fam_after), len(read_bim_ids(f"{stage2_prefix}.bim")),
                         {"mind": args.mind},
                         parse_plink_log_reasons(f"{stage2_prefix}.log"))

    # ---- stage3: 最终位点缺失 / MAF / MAC ----
    stage3_prefix = out_dir / "stage3"
    code = run_plink2(
        ["--bfile", str(stage2_prefix), "--geno", str(args.geno),
         "--maf", str(args.maf), "--mac", str(args.mac),
         *chr_set_args(args.chr_count),
         "--make-bed", "--out", str(stage3_prefix)],
        log_path, args.threads, args.memory_mb,
    )
    if code != 0:
        fail(f"stage3（geno/maf/mac={args.geno}/{args.maf}/{args.mac}）plink2 失败"
             f"（退出码 {code}），日志: {log_path}")
    bim2, bim3 = read_bim_ids(f"{stage2_prefix}.bim"), read_bim_ids(f"{stage3_prefix}.bim")
    removed_v3 = sorted(set(bim2) - set(bim3))
    write_tsv(out_dir / "removed_variants_stage3.tsv",
              ["variant_id"], [[v] for v in removed_v3])
    s3 = stage_stats_row("stage3_geno_maf_mac", s2["n_samples_after"], s2["n_variants_after"],
                         len(read_fam_ids(f"{stage3_prefix}.fam")), len(bim3),
                         {"geno": args.geno, "maf": args.maf, "mac": args.mac},
                         parse_plink_log_reasons(f"{stage3_prefix}.log"))

    # ---- 各阶段独立统计文件（与 final_stats 对账）----
    for stage_idx, stage_row in ((1, s1), (2, s2), (3, s3)):
        write_tsv(out_dir / f"stage{stage_idx}.stats.tsv", STATS_HEADER, [
            [stage_row["stage"], stage_row["n_samples_before"],
             stage_row["n_variants_before"], stage_row["n_samples_removed"],
             stage_row["n_variants_removed"], stage_row["n_samples_after"],
             stage_row["n_variants_after"], json.dumps(stage_row["params"])]
        ])

    # ---- 阻断条件：QC 后样本过少或标记集为空 ----
    if s3["n_samples_after"] < args.min_samples_after_qc or s3["n_variants_after"] == 0:
        funnel = [s1, s2, s3]
        detail = "\n".join(
            f"  {s['stage']}: 样本 {s['n_samples_before']}→{s['n_samples_after']}"
            f"（剔除 {s['n_samples_removed']}），位点 {s['n_variants_before']}"
            f"→{s['n_variants_after']}（剔除 {s['n_variants_removed']}），参数 {s['params']}"
            for s in funnel
        )
        fail(
            f"QC 后样本数 {s3['n_samples_after']}（阈值 {args.min_samples_after_qc}）"
            f"或标记集 {s3['n_variants_after']} 不满足阻断条件。各阶段流失明细:\n{detail}"
        )

    # ---- final_stats: 最终队列重新统计（与逐阶段计数可对账）----
    cohort_prefix = out_dir / "cohort"
    # 以 stage3 结果作为 cohort（复制三件套）
    import shutil

    for ext in (".bed", ".bim", ".fam"):
        shutil.copy2(f"{stage3_prefix}{ext}", f"{cohort_prefix}{ext}")

    freq_code = run_plink2(
        ["--bfile", str(cohort_prefix), "--freq", *chr_set_args(args.chr_count),
         "--out", str(out_dir / "final_stats")],
        log_path, args.threads, args.memory_mb,
    )
    miss_code = run_plink2(
        ["--bfile", str(cohort_prefix), "--missing", *chr_set_args(args.chr_count),
         "--out", str(out_dir / "final_stats")],
        log_path, args.threads, args.memory_mb,
    )
    if freq_code != 0 or miss_code != 0:
        fail(f"final_stats 重新统计失败（freq={freq_code}, missing={miss_code}）")

    final_summary = {
        "n_samples_final": s3["n_samples_after"],
        "n_variants_final": s3["n_variants_after"],
        "recomputed_on_final_cohort": True,
        "freq_table": "final_stats.afreq",
        "missing_tables": ["final_stats.imiss", "final_stats.lmiss"],
    }

    # final_stats.tsv：可直接阅读的阶段对账汇总
    rows = [
        [s["stage"], s["n_samples_before"], s["n_variants_before"],
         s["n_samples_removed"], s["n_variants_removed"],
         s["n_samples_after"], s["n_variants_after"], json.dumps(s["params"])]
        for s in (s1, s2, s3)
    ]
    rows.append(["final_stats", s3["n_samples_after"], s3["n_variants_after"],
                 0, 0, s3["n_samples_after"], s3["n_variants_after"], "recomputed"])
    write_tsv(out_dir / "final_stats.tsv", STATS_HEADER, rows)

    with open(out_dir / "final_stats.json", "w", encoding="utf-8") as handle:
        json.dump({"stages": [s1, s2, s3], "final": final_summary},
                  handle, indent=2, ensure_ascii=False)

    # ---- qc_audit.json（HWE 默认 off 等审计）----
    audit = {
        "schema": "gwasflow.qc_audit.v1",
        "hwe": "off（自交/近交植物默认关闭 HWE 过滤）",
        "relatedness_removal": "off（高亲缘材料保留，不自动删除；仅标记/展示）",
        "chr_set": f"--chr-set {args.chr_count} no-x no-y no-xy no-mt（非人类语义固定）",
        "nonfounders": True,
        "founder_rule": "纳入全部材料的统计与过滤（家系字段不影响样本计数）",
        "stages": [s1, s2, s3],
        "params": {
            "pre_geno": args.pre_geno, "mind": args.mind, "geno": args.geno,
            "maf": args.maf, "mac": args.mac,
            "min_samples_after_qc": args.min_samples_after_qc,
            "chr_count": args.chr_count,
        },
        "notes": "QC 默认阈值为审阅起点（mind=0.10/geno=0.10/maf=0.05/mac=5），非普适植物标准",
    }
    with open(out_dir / "qc_audit.json", "w", encoding="utf-8") as handle:
        json.dump(audit, handle, indent=2, ensure_ascii=False)

    print(
        f"[genotype_qc] QC 完成：样本 {n_samples}→{s3['n_samples_after']}，"
        f"位点 {n_variants}→{s3['n_variants_after']}，产物目录 {out_dir}"
    )
    return EXIT_OK


# ---- structure 模式：LD 剪枝 + 剪枝集 PCA ----

def run_structure(args, out_dir):
    log_path = out_dir / "plink2.log"

    # 1. LD 剪枝（--indep-pairwise 物理窗口 kb；plink2 约束：kb 窗口 step 必须为 1，
    #    故不传 step 用 plink2 默认值，配置 prune_step 仅适用标记数窗口模式）
    prune_prefix = out_dir / "pca_prune"
    code = run_plink2(
        ["--bfile", args.bfile, "--indep-pairwise",
         f"{args.prune_window_kb}kb", str(args.prune_r2),
         *chr_set_args(args.chr_count),
         "--out", str(prune_prefix)],
        log_path, args.threads, args.memory_mb,
    )
    if code != 0:
        fail(f"LD 剪枝失败（退出码 {code}），日志: {log_path}")

    prune_in = prune_prefix.with_suffix(".prune.in")
    prune_out = prune_prefix.with_suffix(".prune.out")
    n_pruned_in = sum(1 for _ in open(prune_in, encoding="utf-8"))
    n_pruned_out = sum(1 for _ in open(prune_out, encoding="utf-8"))
    n_total = len(read_bim_ids(f"{args.bfile}.bim"))

    warnings = []
    if n_pruned_in == 0:
        warnings.append("剪枝集为空（prune.in 无标记），结构分析不可用；请检查剪枝阈值")
    elif n_pruned_in < 50:
        warnings.append(
            f"剪枝集过小（{n_pruned_in} 个标记，诊断性警告；结构分析属诊断性质，继续执行）"
        )

    # 2. 有效矩阵秩校验（请求 PC 数 >= min(N, M) → 报错）
    n_samples = len(read_fam_ids(f"{args.bfile}.fam"))
    effective_rank_cap = min(n_samples - 1, n_pruned_in) if n_pruned_in > 0 else 0
    if args.pca_components >= effective_rank_cap and effective_rank_cap > 0:
        fail(
            f"请求 PC 数 {args.pca_components} 不小于有效矩阵秩上限 "
            f"min(N-1, M_pruned)={effective_rank_cap}（N={n_samples}，"
            f"M_pruned={n_pruned_in}）；结构分析阶段拒绝"
        )
    if n_pruned_in == 0:
        # 剪枝集为空时 PCA 无法执行：输出空表并告警（合法降级）
        write_tsv(out_dir / "pca.eigenvec", ["#FID", "IID", "PC1"], [])
        with open(out_dir / "pca.eigenval", "w", encoding="utf-8") as handle:
            handle.write("")
        with open(out_dir / "structure_audit.json", "w", encoding="utf-8") as handle:
            json.dump({
                "schema": "gwasflow.structure_audit.v1",
                "prune_params": {
                    "window_kb": args.prune_window_kb, "step": 1,
                    "step_note": "plink2 kb 窗口强制 step=1（配置 prune_step 仅适用标记数窗口模式）",
                    "r2": args.prune_r2,
                },
                "n_prune_in": n_pruned_in, "n_prune_out": n_pruned_out,
                "n_total_variants": n_total,
                "pca_components_requested": args.pca_components,
                "warnings": warnings,
                "pca_status": "skipped_empty_prune_set",
                "usage_isolation": "PCA 使用剪枝标记集；关联分析使用完整 cohort 标记集",
            }, handle, indent=2, ensure_ascii=False)
        print(f"[genotype_qc] 剪枝集为空，PCA 跳过（警告已记录）: {out_dir}")
        return EXIT_OK

    # 3. 剪枝集 PCA（scols=fid 显式保留 FID）
    pca_prefix = out_dir / "pca"
    code = run_plink2(
        ["--bfile", args.bfile, "--extract", str(prune_in),
         "--pca", str(args.pca_components), "scols=fid",
         *chr_set_args(args.chr_count),
         "--out", str(pca_prefix)],
        log_path, args.threads, args.memory_mb,
    )
    if code != 0:
        fail(f"PCA 失败（退出码 {code}），日志: {log_path}")

    with open(out_dir / "structure_audit.json", "w", encoding="utf-8") as handle:
        json.dump({
            "schema": "gwasflow.structure_audit.v1",
            "prune_params": {
                "window_kb": args.prune_window_kb, "step": 1,
                "step_note": "plink2 kb 窗口强制 step=1（配置 prune_step 仅适用标记数窗口模式）",
                "r2": args.prune_r2,
            },
            "pca_params": {"components": args.pca_components, "scols": "fid",
                           "marker_set": "LD-pruned"},
            "n_prune_in": n_pruned_in, "n_prune_out": n_pruned_out,
            "n_total_variants": n_total,
            "warnings": warnings,
            "pca_status": "ok",
            "usage_isolation": "PCA 使用剪枝标记集；关联分析使用完整 cohort 标记集",
        }, handle, indent=2, ensure_ascii=False)

    print(
        f"[genotype_qc] 结构分析完成：剪枝集 {n_pruned_in}/{n_total} 标记，"
        f"PCA {args.pca_components} 维，产物目录 {out_dir}"
    )
    return EXIT_OK


def main(argv=None):
    args = parse_args(argv)
    out_dir = ensure_empty_dir(args.out_dir)

    if args.gate and not Path(args.gate).exists():
        fail(f"门控产物不存在: {args.gate}")

    if not Path(f"{args.bfile}.bim").exists() or not Path(f"{args.bfile}.fam").exists():
        fail(f"输入 BED 前缀不存在: {args.bfile}")

    if args.mode == "qc":
        return run_qc(args, out_dir)
    return run_structure(args, out_dir)


if __name__ == "__main__":
    sys.exit(main())
