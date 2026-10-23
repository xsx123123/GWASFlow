#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GWASFlow 基因型规范化（参数手册 §9.2 + §7 命令链）。

VCF/BCF 路线命令链：
  bcftools view（双等位 SNP + FILTER 策略）
  → bcftools norm --check-ref e（REF 对 FASTA 严格校验）
  → 原始身份导出 → bcftools annotate --set-id（稳定 variant_id）
  → bcftools annotate --rename-chrs（engine 数值染色体）
  → bcftools index --csi
  → plink2 --vcf --make-bed --double-id（FID=IID=样本 ID）

BED 路线：三件套 + 已核验 variant_map（BIM A1/A2 不盲猜 REF/ALT）。

红线：
- 保留原始双等位 SNP，排除多等位/INDEL 并逐类计数（禁止默认拆分）
- 未知染色体报错（列出染色体名与映射表路径，禁止统一转 0）
- REF/FASTA 不一致 → 退出 2 指明位点与双方碱基
- DP/GQ 置缺失：显式启用但 FORMAT 字段缺失 → 阻断（禁止静默跳过）
- 剂量输入（representation != hardcall）显式拒绝；倍性 != 2 拒绝

错误码：0=成功；2=格式/配置/输入错误；外部软件非零退出码透传。
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
VALID_BASES = {"A", "C", "G", "T"}


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="normalize_genotype.py",
        description="GWASFlow 基因型规范化：VCF/BCF 或 BED → 引擎 BED + variant_map + 审计。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--genotype-format", required=True, choices=["vcf", "bcf", "bed"])
    parser.add_argument("--input", default=None, help="VCF/BCF 文件路径")
    parser.add_argument("--bfile", default=None, help="BED 前缀（不含 .bed/.bim/.fam）")
    parser.add_argument("--input-variant-map", default=None,
                        help="BED 路线必填 variant_map.tsv（七列契约；BIM A1/A2 无 REF/ALT 断言）")
    parser.add_argument("--gate", required=True, help="01.gate/input_gate.json（门控产物）")
    parser.add_argument("--samples", required=True, help="samples.csv（FAM 对账）")
    parser.add_argument("--genome-fasta", default=None, help="参考 FASTA（VCF 路线必填）")
    parser.add_argument("--fai", default=None, help="FASTA 索引（VCF 路线必填）")
    parser.add_argument("--chromosome-map", default=None,
                        help="染色体重命名映射 TSV（两列无表头：原始名→数值名，VCF 路线必填）")
    parser.add_argument("--ploidy", type=int, default=2, help="倍性（首版仅 2）")
    parser.add_argument("--representation", default="hardcall",
                        choices=["hardcall"], help="基因型表示（剂量输入拒绝）")
    parser.add_argument("--multiallelic-policy", default="exclude",
                        choices=["exclude"], help="多等位策略（split 属后续工作包）")
    parser.add_argument("--filter-policy", default="pass_only",
                        choices=["pass_only", "pass_or_unfiltered"],
                        help="FILTER 策略（pass_or_unfiltered 的 '.' 为未评估语义，"
                             "不称已通过质量过滤）")
    parser.add_argument("--dp-min", type=int, default=None,
                        help="FORMAT/DP 基因型置缺失阈值（显式启用需字段存在）")
    parser.add_argument("--gq-min", type=int, default=None,
                        help="FORMAT/GQ 基因型置缺失阈值（显式启用需字段存在）")
    parser.add_argument("--chr-count", type=int, default=9,
                        help="染色体对数（plink2 --chr-set，1..95）")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--memory-mb", type=int, default=8000)
    parser.add_argument("--out-dir", required=True, help="输出目录（必须为新目录或空目录）")
    return parser.parse_args(argv)


def fail(message, code=EXIT_INPUT_ERROR):
    print(f"[normalize_genotype ERROR] {message}", file=sys.stderr)
    sys.exit(code)


def run_tool(argv, log_path, env_extra=None, cwd=None):
    """执行外部命令数组；stdout/stderr 写入工具日志；退出码透传。"""
    env = os.environ.copy()
    if env_extra:
        env.update({k: str(v) for k, v in env_extra.items()})
    with open(log_path, "a", encoding="utf-8") as handle:
        handle.write(f"\n=== RUN: {' '.join(map(str, argv))} ===\n")
        handle.flush()
        result = subprocess.run(
            list(map(str, argv)), stdout=handle, stderr=subprocess.STDOUT,
            env=env, cwd=cwd, check=False,
        )
    return result.returncode


def ensure_empty_dir(out_dir):
    path = Path(out_dir).resolve()
    if path.exists():
        if not path.is_dir() or any(path.iterdir()):
            fail(f"输出目录非空或非目录（拒绝覆盖既有产物）: {path}")
    else:
        path.mkdir(parents=True, exist_ok=True)
    return path


def load_gate(gate_path):
    p = Path(gate_path)
    if not p.exists():
        fail(f"门控产物不存在: {p}")
    with open(p, encoding="utf-8") as handle:
        return json.load(handle)


def load_chromosome_map(map_path):
    """加载染色体映射表（两列无表头：原始名 → 数值 engine_chrom）。"""
    p = Path(map_path)
    if not p.exists():
        fail(f"染色体重命名映射表不存在: {p}")
    mapping = {}
    with open(p, newline="", encoding="utf-8") as handle:
        for line_no, row in enumerate(csv.reader(handle, delimiter="\t"), start=1):
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 2:
                fail(f"染色体映射表第 {line_no} 行列数不足: {row}")
            orig, engine = row[0].strip(), row[1].strip()
            if orig in mapping:
                fail(f"染色体映射表原始名重复: '{orig}'（第 {line_no} 行）")
            try:
                engine_value = int(engine)
            except ValueError:
                fail(
                    f"染色体映射表第 {line_no} 行 engine 值非整数: '{engine}'"
                )
            mapping[orig] = engine_value
    if not mapping:
        fail(f"染色体映射表为空: {p}")
    return mapping


def load_samples_fid_iid(samples_csv):
    pairs = set()
    with open(samples_csv, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            pairs.add((row.get("fid", "").strip(), row.get("iid", "").strip()))
    return pairs


def classify_variant_records(vcf_path, bcftools, log_path):
    """遍历原始 VCF 记录并分类计数（审计：双等位 SNP/多等位/INDEL/其他）。"""
    argv = [bcftools, "query", "-f", "%CHROM\\t%REF\\t%ALT\\n", str(vcf_path)]
    env = os.environ.copy()
    counts = {
        "total_records": 0,
        "biallelic_snp": 0,
        "multiallelic": 0,
        "indel": 0,
        "other": 0,
    }
    chromosomes = set()
    proc = subprocess.run(
        argv, capture_output=True, text=True, env=env, check=False
    )
    if proc.returncode != 0:
        with open(log_path, "a", encoding="utf-8") as handle:
            handle.write(proc.stderr)
        fail(f"bcftools query 失败（退出码 {proc.returncode}）: {vcf_path}")
    for line in proc.stdout.splitlines():
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            counts["other"] += 1
            continue
        chrom, ref, alt = parts[0], parts[1].upper(), parts[2].upper()
        chromosomes.add(chrom)
        counts["total_records"] += 1
        alts = alt.split(",")
        if len(alts) > 1:
            counts["multiallelic"] += 1
        elif len(ref) == 1 and len(alts[0]) == 1 and ref in VALID_BASES and alts[0] in VALID_BASES:
            counts["biallelic_snp"] += 1
        elif len(ref) != len(alts[0]):
            counts["indel"] += 1
        else:
            counts["other"] += 1
    return counts, chromosomes


def check_format_field(vcf_path, field, bcftools):
    """检测 VCF header 是否声明 FORMAT/<field>。"""
    argv = [bcftools, "view", "-h", str(vcf_path)]
    proc = subprocess.run(argv, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        return None
    return f"FORMAT/{field}" in proc.stdout


def load_fam_ids(fam_path):
    pairs = []
    with open(fam_path, encoding="utf-8") as handle:
        for line in handle:
            parts = line.split()
            if len(parts) >= 2:
                pairs.append((parts[0], parts[1]))
    return pairs


# ---- VCF/BCF 路线 ----

def normalize_vcf(args, out_dir, audit):
    bcftools = "bcftools"
    plink2 = "plink2"
    log_path = out_dir / "normalize.log"

    vcf_in = Path(args.input)
    if not vcf_in.exists():
        fail(f"基因型输入文件不存在: {vcf_in}")

    # 0. 染色体映射校验（未知染色体报错，禁止统一转 0）
    mapping = load_chromosome_map(args.chromosome_map)
    counts, chromosomes = classify_variant_records(vcf_in, bcftools, log_path)
    unknown = sorted(c for c in chromosomes if c not in mapping)
    if unknown:
        fail(
            f"输入包含染色体映射表之外的染色体: {unknown}；"
            f"映射表路径: {args.chromosome_map}（禁止将未知染色体统一转换为 0）"
        )
    audit["variant_filter_counts"] = counts
    audit["chromosome_mapping"] = {c: mapping[c] for c in sorted(chromosomes)}

    # 0.5 DP/GQ 置缺失（先置缺失再进入后续过滤；字段缺失时阻断）
    dp_gq = {
        "dp_min": {"value": args.dp_min, "executed": False, "reason": ""},
        "gq_min": {"value": args.gq_min, "executed": False, "reason": ""},
    }
    current_input = vcf_in
    for field, threshold, key in (("DP", args.dp_min, "dp_min"), ("GQ", args.gq_min, "gq_min")):
        if threshold is None:
            dp_gq[key]["reason"] = "未启用（未提供阈值）"
            continue
        if not check_format_field(current_input, field, bcftools):
            fail(
                f"显式启用 {key}={threshold} 但 FORMAT/{field} 字段不存在"
                f"（禁止静默跳过该过滤）: {current_input}"
            )
        filtered = out_dir / f"dp_gq_filtered_{field.lower()}.vcf.gz"
        argv = [
            bcftools, "filter", "-S", ".",
            "-e", f"FMT/{field}<{threshold} & FMT/{field}!=.",
            "-O", "z", "-o", str(filtered), str(current_input),
        ]
        code = run_tool(argv, log_path)
        if code != 0:
            fail(f"bcftools filter FORMAT/{field} 置缺失失败（退出码 {code}）")
        dp_gq[key]["executed"] = True
        dp_gq[key]["reason"] = f"基因型级置缺失（FMT/{field}<{threshold}）"
        current_input = filtered
    audit["dp_gq_filters"] = dp_gq

    # 1. 双等位 SNP + FILTER 策略
    if args.filter_policy == "pass_only":
        apply_filters = "PASS"
        filter_semantics = "仅保留 FILTER=PASS"
    else:
        apply_filters = "PASS,."
        filter_semantics = "保留 PASS 与未评估 '.'（'.' 表示未评估，不称已通过质量过滤）"
    biallelic = out_dir / "biallelic.vcf.gz"
    argv = [
        bcftools, "view",
        "--min-alleles", "2", "--max-alleles", "2", "--types", "snps",
        "--apply-filters", apply_filters,
        "-O", "z", "-o", str(biallelic), str(current_input),
    ]
    code = run_tool(argv, log_path)
    if code != 0:
        fail(f"bcftools view（双等位 SNP 筛选）失败（退出码 {code}）")
    audit["filter_policy"] = {"policy": args.filter_policy, "semantics": filter_semantics}

    # 2. REF 对 FASTA 严格校验（--check-ref e）
    normalized = out_dir / "normalized.vcf.gz"
    argv = [
        bcftools, "norm",
        "--fasta-ref", args.genome_fasta, "--check-ref", "e",
        "-O", "z", "-o", str(normalized), str(biallelic),
    ]
    code = run_tool(argv, log_path)
    if code != 0:
        mismatch_lines = []
        with open(log_path, encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if "REF_MISMATCH" in line or "ref" in line.lower() and "mismatch" in line.lower():
                    mismatch_lines.append(line.strip())
        detail = "; ".join(mismatch_lines[-5:]) if mismatch_lines else "详见工具日志"
        fail(
            f"REF 与 FASTA 不一致（--check-ref e 拒绝）：{detail}；"
            f"完整日志: {log_path}"
        )

    # 3. 原始身份导出（variant_map 数据源）
    original_rows = []
    argv = [bcftools, "query", "-f", "%ID\\t%CHROM\\t%POS\\t%REF\\t%ALT\\n", str(normalized)]
    proc = subprocess.run(argv, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        fail("bcftools query（原始身份导出）失败")
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) == 5:
            original_rows.append(parts)

    # 4. 稳定 variant_id 置换
    identified = out_dir / "identified.vcf.gz"
    argv = [
        bcftools, "annotate",
        "--set-id", "%CHROM:%POS:%REF:%FIRST_ALT",
        "-O", "z", "-o", str(identified), str(normalized),
    ]
    code = run_tool(argv, log_path)
    if code != 0:
        fail(f"bcftools annotate --set-id 失败（退出码 {code}）")

    # 5. 染色体重命名（engine 数值编码）
    engine = out_dir / "engine.vcf.gz"
    argv = [
        bcftools, "annotate",
        "--rename-chrs", args.chromosome_map,
        "-O", "z", "-o", str(engine), str(identified),
    ]
    code = run_tool(argv, log_path)
    if code != 0:
        fail(f"bcftools annotate --rename-chrs 失败（退出码 {code}）")

    # 6. 索引
    code = run_tool([bcftools, "index", "--csi", str(engine)], log_path)
    if code != 0:
        fail(f"bcftools index 失败（退出码 {code}）")

    # 7. plink2 导入（--double-id FID=IID、半分型报错、非人类染色体集）
    raw_prefix = out_dir / "raw_engine"
    chr_set_args = [
        "--chr-set", str(args.chr_count), "no-x", "no-y", "no-xy", "no-mt",
    ]
    argv = [
        plink2, "--vcf", str(engine),
        "--double-id", "--vcf-half-call", "error",
        *chr_set_args,
        "--threads", str(args.threads),
        "--memory", str(args.memory_mb),
        "--make-bed", "--out", str(raw_prefix),
    ]
    code = run_tool(argv, log_path)
    if code != 0:
        fail(f"plink2 --vcf 导入失败（退出码 {code}）；半分型/倍性问题详见日志 {log_path}")

    # 8. FAM 与 samples.csv 对账
    sample_pairs = load_samples_fid_iid(args.samples)
    fam_pairs = load_fam_ids(raw_prefix.with_suffix(".fam"))
    fam_set = set(fam_pairs)
    unknown_in_fam = fam_set - sample_pairs
    if unknown_in_fam:
        fail(
            f"导入 FAM 中存在 samples.csv 未登记的身份: "
            f"{sorted(unknown_in_fam)}（--double-id 要求 FID=IID=VCF 样本 ID）"
        )
    missing_in_fam = sample_pairs - fam_set
    if missing_in_fam:
        audit["samples_crosscheck"] = {
            "fam_n": len(fam_set),
            "samples_csv_n": len(sample_pairs),
            "registered_not_in_genotype": sorted(missing_in_fam),
            "note": "samples.csv 登记但基因型中不存在的材料（后续表型准备按 QC 队列处理）",
        }
    else:
        audit["samples_crosscheck"] = {"fam_n": len(fam_set), "samples_csv_n": len(sample_pairs)}

    # 9. variant_map.tsv 生成（七列带表头）
    variant_map_path = out_dir / "variant_map.tsv"
    with open(variant_map_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "variant_id", "original_id", "chrom", "engine_chrom", "pos", "ref", "alt",
        ])
        for original_id, chrom, pos, ref, alt in original_rows:
            variant_id = f"{chrom}:{pos}:{ref}:{alt}"
            writer.writerow([
                variant_id, original_id if original_id != "." else ".",
                chrom, mapping[chrom], pos, ref, alt,
            ])

    # 10. 引擎产物就位（raw_engine.* 已在 out_dir）
    return {
        "bed": str(raw_prefix.with_suffix(".bed")),
        "bim": str(raw_prefix.with_suffix(".bim")),
        "fam": str(raw_prefix.with_suffix(".fam")),
        "variant_map": str(variant_map_path),
    }


# ---- BED 路线 ----

def normalize_bed(args, out_dir, audit):
    bfile = Path(args.bfile)
    for ext in (".bed", ".bim", ".fam"):
        if not Path(str(bfile) + ext).exists():
            fail(f"BED 三件套不齐备，缺少 {str(bfile) + ext}")

    variant_map_in = Path(args.input_variant_map) if args.input_variant_map else None
    if not variant_map_in or not variant_map_in.exists():
        fail(
            "BED 入口必须伴随已核验的 variant_map（--input-variant-map）；"
            "BIM 的 A1/A2 不能独立证明 REF/ALT"
        )

    # variant_map 与 BIM 对账（身份/坐标/等位基因）
    bim_records = {}
    with open(str(bfile) + ".bim", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            parts = line.split()
            if len(parts) < 6:
                fail(f"BIM 第 {line_no} 行列数不足")
            chrom, vid, cm, pos, a1, a2 = parts[0], parts[1], parts[2], parts[3], parts[4], parts[5]
            bim_records[vid] = (chrom, pos, a1, a2)

    mapped = {}
    with open(variant_map_in, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for col in ("variant_id", "original_id", "chrom", "engine_chrom", "pos", "ref", "alt"):
            if col not in (reader.fieldnames or []):
                fail(f"variant_map 缺失列: {col}")
        for line_no, row in enumerate(reader, start=2):
            vid = row["variant_id"]
            if vid in mapped:
                fail(f"variant_map variant_id 重复: '{vid}'（第 {line_no} 行）")
            mapped[vid] = row

    missing_in_map = set(bim_records) - set(mapped)
    if missing_in_map:
        fail(
            f"BIM 中存在 variant_map 未覆盖的位点 ID（前 10 个）: "
            f"{sorted(missing_in_map)[:10]}"
        )
    for vid, (chrom, pos, a1, a2) in bim_records.items():
        row = mapped[vid]
        if str(row["engine_chrom"]) != str(chrom) or str(row["pos"]) != str(pos):
            fail(
                f"variant_map 与 BIM 坐标不一致: {vid} "
                f"(map: {row['engine_chrom']}:{row['pos']}, bim: {chrom}:{pos})"
            )
        if {row["ref"], row["alt"]} != {a1, a2}:
            fail(
                f"variant_map 与 BIM 等位基因不一致: {vid} "
                f"(map: {row['ref']}/{row['alt']}, bim: {a1}/{a2})"
            )

    # 复制三件套为 raw_engine.*
    import shutil

    raw_prefix = out_dir / "raw_engine"
    for ext in (".bed", ".bim", ".fam"):
        shutil.copy2(str(bfile) + ext, str(raw_prefix) + ext)

    # FAM 与 samples.csv 对账
    sample_pairs = load_samples_fid_iid(args.samples)
    fam_pairs = set(load_fam_ids(raw_prefix.with_suffix(".fam")))
    unknown_in_fam = fam_pairs - sample_pairs
    if unknown_in_fam:
        fail(f"导入 FAM 中存在 samples.csv 未登记的身份: {sorted(unknown_in_fam)}")
    audit["samples_crosscheck"] = {"fam_n": len(fam_pairs), "samples_csv_n": len(sample_pairs)}

    # variant_map 复制为输出产物
    variant_map_out = out_dir / "variant_map.tsv"
    shutil.copy2(variant_map_in, variant_map_out)

    audit["bed_route"] = {
        "source_bfile": str(bfile),
        "variant_map_source": str(variant_map_in),
        "n_bim_variants": len(bim_records),
        "identity_check": "BIM/variant_map 坐标与等位基因全量对账通过",
    }
    return {
        "bed": str(raw_prefix) + ".bed",
        "bim": str(raw_prefix) + ".bim",
        "fam": str(raw_prefix) + ".fam",
        "variant_map": str(variant_map_out),
    }


def main(argv=None):
    args = parse_args(argv)
    out_dir = ensure_empty_dir(args.out_dir)
    load_gate(args.gate)

    if args.ploidy != 2:
        fail(f"首版仅支持二倍体（ploidy=2），当前: {args.ploidy}")
    # representation 参数 choices 已限定 hardcall（剂量输入拒绝）

    audit = {
        "schema": "gwasflow.normalize_audit.v1",
        "genotype_format": args.genotype_format,
        "ploidy": args.ploidy,
        "representation": args.representation,
        "multiallelic_policy": args.multiallelic_policy,
        "chr_count": args.chr_count,
    }

    if args.genotype_format in ("vcf", "bcf"):
        if not args.input or not args.genome_fasta or not args.chromosome_map:
            fail("VCF/BCF 路线必须提供 --input/--genome-fasta/--chromosome-map")
        products = normalize_vcf(args, out_dir, audit)
    else:
        if not args.bfile:
            fail("BED 路线必须提供 --bfile")
        products = normalize_bed(args, out_dir, audit)

    audit["products"] = products
    with open(out_dir / "normalize_audit.json", "w", encoding="utf-8") as handle:
        json.dump(audit, handle, indent=2, ensure_ascii=False)

    print(f"[normalize_genotype] 规范化完成，产物目录: {out_dir}")
    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
