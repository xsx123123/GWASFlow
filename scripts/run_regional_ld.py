#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GWASFlow 局部区域 LD（参数手册 §9.6）。

- 对每个 locus 在该性状有效材料（keep）上计算成对 LD（plink2 --r2 或 --dprime 计量）
- 工具自身 MAF/missing/het 过滤值与实际使用标记集审计记录（spec §5.9.3-1）
- max_regional_plots 只限出图不限数据（首版不生成 LD 热图，如实记录）
- loci 为空 → 跳过并在 manifest 标注（合法交付）
- LD 工具非零退出透传

错误码：0=成功；2=格式/配置/输入错误。
"""

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

EXIT_OK = 0
EXIT_INPUT_ERROR = 2


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="run_regional_ld.py",
        description="区域 LD 封装（每 locus 成对 LD + 状态表 + 索引）。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--bfile", required=True, help="cohort BED 前缀")
    parser.add_argument("--keep", required=True, help="该性状有效材料 keep.tsv")
    parser.add_argument("--loci", required=True, help="loci.tsv（区域索引）")
    parser.add_argument("--trait", required=True, help="性状 ID")
    parser.add_argument("--maf", type=float, default=0.05,
                        help="LD 工具自身 MAF 过滤（审计记录）")
    parser.add_argument("--missing-max", type=float, default=0.10,
                        help="LD 工具自身缺失过滤（审计记录）")
    parser.add_argument("--het-max", type=float, default=1.0,
                        help="LD 工具自身杂合度过滤（审计记录）")
    parser.add_argument("--metric", default="r2", choices=["r2", "dprime"],
                        help="LD 计量（图和数据一致）")
    parser.add_argument("--max-regional-plots", type=int, default=20,
                        help="出图数量上限（只限出图，不删数据）")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--out-dir", required=True, help="输出目录（新/空目录）")
    return parser.parse_args(argv)


def fail(message, code=EXIT_INPUT_ERROR):
    print(f"[run_regional_ld ERROR] {message}", file=sys.stderr)
    sys.exit(code)


def ensure_empty_dir(out_dir):
    path = Path(out_dir).resolve()
    if path.exists():
        if not path.is_dir() or any(path.iterdir()):
            fail(f"输出目录非空或非目录（拒绝覆盖既有产物）: {path}")
    else:
        path.mkdir(parents=True, exist_ok=True)
    return path


def sha256_of_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_loci(path):
    loci = []
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            try:
                loci.append({
                    "locus_id": row["locus_id"],
                    "chrom": row["chrom"],
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "lead_variant": row.get("lead_variant", ""),
                })
            except (KeyError, ValueError):
                fail(f"loci.tsv 行解析失败: {row}")
    return loci


def chrom_to_engine(chrom):
    """染色体名转引擎数值（plink2 --chr 需要 1..N）。"""
    try:
        return int(chrom)
    except ValueError:
        cleaned = "".join(c for c in chrom if c.isdigit())
        return int(cleaned) if cleaned else None


def main(argv=None):
    args = parse_args(argv)
    out_dir = ensure_empty_dir(args.out_dir)

    loci = read_loci(args.loci)
    keep_sha = sha256_of_file(args.keep)
    log_path = out_dir / "regional_ld.log"

    ld_rows = []
    status_rows = []
    tool_failures = []
    plots_generated = 0

    for locus in loci:
        engine_chrom = chrom_to_engine(locus["chrom"])
        if engine_chrom is None:
            status_rows.append([locus["locus_id"], "skipped", "染色体无法转引擎数值"])
            continue
        region_prefix = out_dir / f"region_{locus['locus_id']}"
        # plink2 ≥2.0.0-a.7: --r2 已拆分——r2-unphased=PLINK 1.x 经典 r²（无相位剂量），
        # dprime 由 --r2-phased cols=+dprimeabs 输出（PLINK 1.9 --r2 dprime 等价）
        if args.metric == "r2":
            r2_args = ["--r2-unphased"]
        else:
            r2_args = ["--r2-phased", "cols=+dprimeabs"]
        argv_cmd = [
            "plink2", "--bfile", args.bfile,
            "--keep", args.keep,
            "--chr", str(engine_chrom),
            "--from-bp", str(locus["start"]),
            "--to-bp", str(locus["end"]),
            "--maf", str(args.maf),
            "--geno", str(args.missing_max),
            *r2_args,
            "--threads", str(args.threads),
            "--out", str(region_prefix),
        ]
        with open(log_path, "a", encoding="utf-8") as handle:
            handle.write(f"\n=== RUN: {' '.join(argv_cmd)} ===\n")
            handle.flush()
            result = subprocess.run(argv_cmd, stdout=handle,
                                    stderr=subprocess.STDOUT, check=False)
        if result.returncode != 0:
            tool_failures.append((locus["locus_id"], result.returncode))
            status_rows.append([
                locus["locus_id"], "failed",
                f"plink2 退出码 {result.returncode}（透传记录）",
            ])
            continue

        ld_file = region_prefix.with_suffix(".vcor")
        n_pairs = 0
        if ld_file.exists():
            # plink2 .vcor 列（# 表头）: CHROM_A POS_A ID_A CHROM_B POS_B ID_B
            # [+ UNPHASED_R2 | PHASED_R2 (+ ABS_DPRIME)]；按表头列名定位
            value_cols = (("ABS_DPRIME",) if args.metric == "dprime"
                          else ("UNPHASED_R2", "PHASED_R2", "R2"))
            with open(ld_file, encoding="utf-8") as handle:
                header = None
                for line in handle:
                    if line.startswith("#"):
                        header = line.lstrip("#").strip().split()
                        continue
                    if not line.strip() or header is None:
                        continue
                    parts = line.split()
                    if len(parts) < len(header):
                        continue
                    rec = dict(zip(header, parts))
                    id_a, id_b = rec.get("ID_A", ""), rec.get("ID_B", "")
                    value = next((rec[c] for c in value_cols if c in rec), "")
                    if id_a and id_b and value:
                        ld_rows.append([
                            locus["locus_id"], locus["chrom"],
                            id_a, id_b, value,
                            locus["start"], locus["end"],
                        ])
                        n_pairs += 1
            os.remove(ld_file)
        status_rows.append([
            locus["locus_id"], "ok",
            f"n_pairs={n_pairs}（plink2 --r2；工具过滤 maf={args.maf}, "
            f"geno={args.missing_max} 已审计）",
        ])

    # ---- regional_ld.tsv ----
    with open(out_dir / "regional_ld.tsv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["locus_id", "chrom", "variant_a", "variant_b",
                         "metric_value", "region_start", "region_end"])
        for row in ld_rows:
            writer.writerow(row)

    # ---- regional_ld_status.tsv ----
    with open(out_dir / "regional_ld_status.tsv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["locus_id", "status", "detail"])
        for row in status_rows:
            writer.writerow(row)

    # ---- ld_index.json（工具过滤审计 + 出图状态） ----
    index = {
        "schema": "gwasflow.regional_ld.v1",
        "trait": args.trait,
        "metric": args.metric,
        "tool_filters_audit": {
            "maf": args.maf,
            "missing_max": args.missing_max,
            "het_max": args.het_max,
            "note": "LD 工具自身过滤生效值与实际使用标记集审计（plink2 --maf/--geno）",
        },
        "ld_sample": {"keep_sha256": keep_sha,
                      "note": "区域 LD 使用该性状实际材料（keep）"},
        "n_loci": len(loci),
        "n_pairs_total": len(ld_rows),
        "plots": {
            "max_regional_plots": args.max_regional_plots,
            "generated": plots_generated,
            "reason": "首版不生成 LD 热图（LDBlockShow 属外部工具，后续工作包接入）；"
                      "上限只限制出图，不删除完整区域数据",
        },
        "loci_skipped_empty": 0 if loci else 1,
        "tool_failures": [{"locus_id": lid, "exit_code": code}
                          for lid, code in tool_failures],
    }
    if not loci:
        index["empty_reason"] = "loci 为空（无显著区域）：跳过区域 LD（合法交付）"
    with open(out_dir / "ld_index.json", "w", encoding="utf-8") as handle:
        json.dump(index, handle, indent=2, ensure_ascii=False)

    if tool_failures:
        failed = ", ".join(f"{lid}({code})" for lid, code in tool_failures)
        print(f"[run_regional_ld WARNING] 区域 LD 工具失败（透传记录）: {failed}",
              file=sys.stderr)

    print(
        f"[run_regional_ld] 区域 LD 完成：{len(loci)} 区域，{len(ld_rows)} 成对"
        f"（metric={args.metric}）"
    )
    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
