#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GWASFlow 统一关联统计（参数手册 §5 契约重实现 + 精确 MAC 注入）。

五件套输出：gwas.tsv / significant.tsv / invalid_variants.tsv /
unreported_variants.tsv / metadata.json。

红线：
- 身份强校验：assoc 的 chr/ps 必须对应 variant_map 的 engine_chrom/pos，等位基因集合
  必须匹配；结构错误直接退出 2（禁止静默跳过）
- 方向转换四同步：beta×-1、eaf=1-eaf、等位基因对调、SE 与 P 不变（禁止从 P 推造效应）
- Decimal 极小 P（1e-1000 原值保留）；P=0 → p_underflow、neg_log10_p=NA
- Bonferroni 阈值 = alpha/(family_size×M)，M=该性状有效 P 个数（记录实际检验数）
- --trait-mac 注入精确 MAC（消除 af×N 舍入值）
- 状态五级：prepared|software_exited_zero|parsed|partial|complete

错误码：0=成功；2=格式/配置/输入错误；3=已导出结果但含数值异常
（--allow-partial 仅影响退出码，状态仍 partial）。
"""

import argparse
import csv
import hashlib
import json
import sys
from decimal import Decimal, InvalidOperation
from pathlib import Path

EXIT_OK = 0
EXIT_INPUT_ERROR = 2
EXIT_PARTIAL = 3

GWAS_COLUMNS = [
    "schema_version", "analysis_id", "trait", "model", "test", "genome_version",
    "chrom", "pos", "variant_id", "original_id", "ref", "alt",
    "effect_allele", "other_allele", "beta", "se", "effect_scale",
    "p_value", "p_value_text", "neg_log10_p", "eaf", "maf", "mac", "n_obs",
    "status", "filter_reason", "significant",
]

SCHEMA_VERSION = "gwasflow.gwas.v1"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="harmonize_gemma.py",
        description="核验运行/性状/assoc 身份、坐标和等位基因，导出统一统计与 Bonferroni 结果。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--assoc", required=True, help="GEMMA 原始 assoc 表")
    parser.add_argument("--variant-map", required=True,
                        help="variant_map.tsv（variant_id,chrom,engine_chrom,pos,ref,alt）")
    parser.add_argument("--sample-metadata", required=True, help="trait 准备 metadata.json")
    parser.add_argument("--run-metadata", required=True, help="GEMMA run_metadata.json（防错配）")
    parser.add_argument("--trait-mac", required=True, help="精确 MAC 表（trait_mac.tsv）")
    parser.add_argument("--analysis-id", required=True, help="一次方法/样本集合确定的分析 ID")
    parser.add_argument("--genome-version", required=True, help="当前参考组装标识")
    parser.add_argument("--effect-scale", required=True, help="效应单位（如 cm、log(cm)；不由脚本猜测）")
    parser.add_argument("--trait", required=True, help="性状 ID")
    parser.add_argument("--p-column", default="p_wald",
                        help="主检验 P 值列（p_wald/p_lrt/p_score；必须属于运行记录的检验）")
    parser.add_argument("--target-effect-allele", default="reported",
                        choices=["reported", "alt"], help="保留 allele1 方向或转换到 ALT")
    parser.add_argument("--alpha", default="0.05", help="(0,1) Decimal 字符串")
    parser.add_argument("--family-size", type=int, default=1,
                        help="预声明主性状数（阈值 alpha/(family_size×M)；部分性状失败不得事后缩小）")
    parser.add_argument("--allow-partial", action="store_true",
                        help="数值异常仍标 partial，但允许返回 0")
    parser.add_argument("--out-dir", required=True, help="标准化结果目录（新/空目录）")
    return parser.parse_args(argv)


def fail(message, code=EXIT_INPUT_ERROR):
    print(f"[harmonize_gemma ERROR] {message}", file=sys.stderr)
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


def load_json(path, label):
    p = Path(path)
    if not p.exists():
        fail(f"{label} 不存在: {p}")
    with open(p, encoding="utf-8") as handle:
        return json.load(handle)


def parse_decimal(text):
    try:
        return Decimal(str(text).strip())
    except InvalidOperation:
        return None


def neg_log10_decimal(value: Decimal):
    """Decimal 精确 -log10(p)；p=0 返回 None（禁止编造）。"""
    if value == 0:
        return None
    return -value.log10()


def format_decimal(value):
    if value is None:
        return "NA"
    return format(value, "g") if value == value else "NA"




def _export_failed_empty_set(args, out_dir, run_meta, sample_meta):
    """上游软件失败的性状：合法空表五件套 + 诚实状态（partial 不写 complete）。"""
    run_status = run_meta.get("status", "failed_software")
    with open(out_dir / "gwas.tsv", "w", newline="", encoding="utf-8") as handle:
        csv.writer(handle, delimiter="\t").writerow(GWAS_COLUMNS)
    with open(out_dir / "significant.tsv", "w", newline="", encoding="utf-8") as handle:
        csv.writer(handle, delimiter="\t").writerow(GWAS_COLUMNS)
    with open(out_dir / "invalid_variants.tsv", "w", newline="", encoding="utf-8") as handle:
        csv.writer(handle, delimiter="\t").writerow(
            ["variant_id", "assoc_line", "status", "filter_reason"])
    with open(out_dir / "unreported_variants.tsv", "w", newline="", encoding="utf-8") as handle:
        csv.writer(handle, delimiter="\t").writerow(
            ["variant_id", "chrom", "pos", "reason_note"])
    metadata = {
        "schema": "gwasflow.harmonized.v1",
        "analysis_id": args.analysis_id,
        "trait": args.trait,
        "test": run_meta.get("test", "unknown"),
        "p_column": args.p_column,
        "target_effect_allele": args.target_effect_allele,
        "status": run_status,
        "counts": {
            "n_reported": 0, "n_valid_p": 0, "n_invalid": 0,
            "n_significant": 0, "n_unreported": 0, "n_total_variant_map": 0,
        },
        "multiple_testing": {
            "method": "Bonferroni", "alpha": str(args.alpha),
            "family_size": args.family_size, "m_tests": 0, "threshold": "NA",
            "note": "上游软件未成功退出，无有效检验；M=0 不参与阈值计算",
        },
        "effect_scale": args.effect_scale,
        "genome_version": args.genome_version,
        "failure": {
            "upstream_status": run_status,
            "upstream_exit_code": run_meta.get("exit_code"),
            "note": "GEMMA 非零退出；整体状态 partial（其余性状不受影响）",
        },
        "fingerprint": {
            "run_metadata_sha256": sha256_of_file(args.run_metadata),
        },
        "run_status_binding": run_status,
        "allow_partial_flag": bool(args.allow_partial),
        "schema_version_field": SCHEMA_VERSION,
    }
    with open(out_dir / "metadata.json", "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, ensure_ascii=False)


def main(argv=None):
    args = parse_args(argv)
    out_dir = ensure_empty_dir(args.out_dir)

    # ---- 前置校验：run-metadata 与 sample-metadata 同一性状/同一次运行 ----
    run_meta = load_json(args.run_metadata, "run_metadata")
    sample_meta = load_json(args.sample_metadata, "sample_metadata")

    if run_meta.get("trait") != args.trait or sample_meta.get("trait") != args.trait:
        fail(
            f"元数据性状不匹配：run_metadata.trait={run_meta.get('trait')!r}, "
            f"sample_metadata.trait={sample_meta.get('trait')!r}, "
            f"参数 --trait={args.trait!r}"
        )
    bound_trait_meta = (run_meta.get("bindings", {}) or {}).get("trait_metadata", {})
    actual_sha = sha256_of_file(args.sample_metadata)
    if bound_trait_meta.get("sha256") not in (None, actual_sha):
        fail(
            f"run_metadata 绑定的 trait metadata 指纹 "
            f"{bound_trait_meta.get('sha256')[:12]}… 与实际文件 {actual_sha[:12]}… "
            f"不一致（防错配校验失败）"
        )
    run_test = run_meta.get("test", "wald")
    p_column_by_test = {"wald": "p_wald", "lrt": "p_lrt", "score": "p_score"}
    allowed_p_columns = (
        ["p_wald", "p_lrt", "p_score"] if run_test == "all"
        else [p_column_by_test.get(run_test, "p_wald")]
    )
    if args.p_column not in allowed_p_columns:
        fail(
            f"--p-column {args.p_column} 不属于运行记录的检验（run test={run_test}，"
            f"允许: {allowed_p_columns}）"
        )
    run_failed = run_meta.get("status") not in ("software_exited_zero",)
    if run_failed:
        # 运行失败传导：输出合法空表 + 诚实状态（场景 9；数值失败不伪装完整成功）
        _export_failed_empty_set(args, out_dir, run_meta, sample_meta)
        print(
            f"[harmonize_gemma] 性状 '{args.trait}' 上游软件未成功退出"
            f"（status={run_meta.get('status')!r}, exit_code={run_meta.get('exit_code')}）；"
            f"已导出合法空表与 failed 状态（failure 传导，供 project_summary 记录）"
        )
        return EXIT_OK

    alpha = parse_decimal(args.alpha)
    if alpha is None or not (Decimal(0) < alpha < Decimal(1)):
        fail(f"--alpha 必须在 (0,1)：{args.alpha}")
    if args.family_size < 1:
        fail(f"--family-size 必须为正整数：{args.family_size}")

    # ---- variant_map 加载（身份对账基准：规范化时点原始身份）----
    variant_map = {}
    engine_index = {}
    with open(args.variant_map, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for col in ("variant_id", "original_id", "chrom", "engine_chrom", "pos", "ref", "alt"):
            if col not in (reader.fieldnames or []):
                fail(f"variant_map 缺失列: {col}")
        for line_no, row in enumerate(reader, start=2):
            vid = row["variant_id"].strip()
            if vid in variant_map:
                fail(f"variant_map variant_id 重复: '{vid}'（第 {line_no} 行）")
            entry = {
                "variant_id": vid,
                "original_id": row["original_id"].strip(),
                "chrom": row["chrom"].strip(),
                "engine_chrom": row["engine_chrom"].strip(),
                "pos": int(row["pos"]),
                "ref": row["ref"].strip().upper(),
                "alt": row["alt"].strip().upper(),
            }
            variant_map[vid] = entry
            engine_index[(entry["engine_chrom"], entry["pos"])] = entry

    # ---- 精确 MAC 加载 ----
    trait_mac = {}
    with open(args.trait_mac, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            vid = row.get("variant_id", "").strip()
            mac_text = (row.get("mac") or "").strip()
            try:
                trait_mac[vid] = int(mac_text)
            except ValueError:
                trait_mac[vid] = None

    # ---- assoc 解析与身份强校验 ----
    assoc_rows = []
    with open(args.assoc, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t", skipinitialspace=True)
        fieldnames = [c.strip() for c in (reader.fieldnames or [])]
        if not fieldnames:
            fail("assoc 文件无表头")
        reader.fieldnames = fieldnames
        for line_no, row in enumerate(reader, start=2):
            row = {k.strip(): (v.strip() if v else "") for k, v in row.items()}
            row["_line_no"] = line_no
            assoc_rows.append(row)

    structural_errors = []
    seen_ids = {}
    harmonized = []
    invalid_rows = []

    for row in assoc_rows:
        line_no = row["_line_no"]
        rs = row.get("rs") or row.get("SNP") or ""
        chrom_raw = row.get("chr") or row.get("CHR") or ""
        a1 = (row.get("allele1") or row.get("A1") or "").upper()
        a2 = (row.get("allele2") or row.get("A2") or "").upper()

        entry = variant_map.get(rs)
        if entry is None:
            structural_errors.append(
                f"assoc 第 {line_no} 行位点 ID '{rs}' 不在 variant_map 中（缺失映射）"
            )
            continue
        if rs in seen_ids:
            structural_errors.append(
                f"assoc 第 {line_no} 行位点 ID 重复: '{rs}'（第 {seen_ids[rs]} 行已出现）"
            )
            continue
        seen_ids[rs] = line_no

        if str(entry["engine_chrom"]) != str(chrom_raw):
            structural_errors.append(
                f"assoc 第 {line_no} 行 '{rs}' 染色体不匹配: assoc chr={chrom_raw!r} vs "
                f"variant_map engine_chrom={entry['engine_chrom']!r}"
            )
            continue
        if {a1, a2} != {entry["ref"], entry["alt"]}:
            structural_errors.append(
                f"assoc 第 {line_no} 行 '{rs}' 等位基因不匹配: assoc {a1}/{a2} vs "
                f"variant_map {entry['ref']}/{entry['alt']}"
            )
            continue

        # ---- 数值解析（beta/se/af/n_obs/P）----
        beta = se = eaf = None
        try:
            beta = float(row.get("beta") or "nan")
            se = float(row.get("se") or "nan")
        except ValueError:
            pass
        try:
            eaf = float(row.get("af") or "nan")
        except ValueError:
            pass
        n_obs = None
        n_miss_text = row.get("n_miss")
        n_obs_text = row.get("n_obs")
        if n_obs_text not in (None, "", "."):
            try:
                n_obs = int(n_obs_text)
            except ValueError:
                pass
        elif n_miss_text not in (None, "", "."):
            try:
                n_analysis = sample_meta.get("n_valid")
                if n_analysis:
                    n_obs = int(n_analysis) - int(n_miss_text)
            except (ValueError, TypeError):
                pass

        p_text = (row.get(args.p_column) or "").strip()
        p_value = parse_decimal(p_text)
        if p_value is None and p_text not in ("", "NA", "."):
            p_value = None

        numeric_ok = (
            beta is not None and beta == beta
            and se is not None and se == se
            and p_value is not None
        )

        # ---- 方向转换（target=alt 且 allele1≠alt：beta×-1、eaf 互补、等位基因对调）----
        effect_allele, other_allele = a1, a2
        final_beta, final_eaf = beta, eaf
        converted = False
        if args.target_effect_allele == "alt" and a1 != entry["alt"]:
            effect_allele, other_allele = a2, a1
            final_beta = -beta if (beta is not None and beta == beta) else beta
            final_eaf = (1.0 - eaf) if (eaf is not None and eaf == eaf) else eaf
            converted = True

        # ---- 极小 P / 下溢保护 ----
        status = "ok"
        filter_reason = ""
        neg_log10 = None
        if numeric_ok and p_value == 0:
            status = "p_underflow"
            filter_reason = "software reported P=0（下溢；原值保留，neg_log10_p=NA）"
            neg_log10 = None
        elif numeric_ok:
            neg_log10 = neg_log10_decimal(p_value)
        if not numeric_ok:
            status = "invalid"
            filter_reason = (
                f"数值解析失败（beta={row.get('beta')!r}, se={row.get('se')!r}, "
                f"p={p_text!r}）"
            )

        maf = None
        if final_eaf is not None and final_eaf == final_eaf:
            maf = min(final_eaf, 1.0 - final_eaf)

        harmonized.append({
            "chrom": entry["chrom"],
            "pos": entry["pos"],
            "variant_id": rs,
            "original_id": entry["original_id"],
            "ref": entry["ref"],
            "alt": entry["alt"],
            "effect_allele": effect_allele,
            "other_allele": other_allele,
            "beta": final_beta,
            "se": se,
            "p_value": p_value,
            "p_value_text": p_text,
            "neg_log10_p": neg_log10,
            "eaf": final_eaf,
            "maf": maf,
            "mac": trait_mac.get(rs),
            "n_obs": n_obs,
            "status": status,
            "filter_reason": filter_reason,
            "converted": converted,
            "line_no": line_no,
        })

    if structural_errors:
        detail = "\n".join(f"  - {e}" for e in structural_errors[:20])
        if len(structural_errors) > 20:
            detail += f"\n  ...（共 {len(structural_errors)} 项冲突）"
        fail(f"assoc 与 variant_map 身份强校验失败（结构错误，禁止静默跳过）:\n{detail}")

    # ---- 排序（chrom, pos, variant_id；小/中规模内存排序路线）----
    def sort_key(rec):
        try:
            chrom_key = int(rec["chrom"])
        except ValueError:
            chrom_key = 10**9
        return (chrom_key, rec["pos"], rec["variant_id"])

    harmonized.sort(key=sort_key)

    # ---- Bonferroni 阈值（M=有效 P 个数；记录实际检验数）----
    valid_p_rows = [r for r in harmonized if r["status"] in ("ok", "p_underflow")]
    m_tests = len(valid_p_rows)
    threshold = (
        alpha / (Decimal(args.family_size) * Decimal(m_tests))
        if m_tests > 0 else None
    )

    invalid_rows = [r for r in harmonized if r["status"] == "invalid"]
    n_invalid = len(invalid_rows)

    # ---- unreported：variant_map 中未出现在 assoc 的 ID ----
    unreported = [
        entry for vid, entry in variant_map.items() if vid not in seen_ids
    ]

    # ---- 状态分级 ----
    run_status = run_meta.get("status", "software_exited_zero")
    if run_status == "software_exited_zero":
        if n_invalid == 0:
            status = "parsed"
        else:
            status = "partial"
    else:
        status = run_status

    # ---- gwas.tsv ----
    def fmt_num(value):
        if value is None:
            return "NA"
        if isinstance(value, float):
            return repr(value) if value == value else "NA"
        if isinstance(value, Decimal):
            return format_decimal(value)
        return str(value)

    significant_rows = []
    with open(out_dir / "gwas.tsv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(GWAS_COLUMNS)
        for r in harmonized:
            significant = False
            if threshold is not None and r["status"] in ("ok", "p_underflow"):
                if r["p_value"] is not None and r["p_value"] != 0:
                    significant = r["p_value"] <= threshold
                elif r["status"] == "p_underflow":
                    significant = True   # P=0 下溢必然小于任何正阈值
            if significant:
                significant_rows.append(r)
            writer.writerow([
                SCHEMA_VERSION, args.analysis_id, args.trait,
                sample_meta.get("primary_model", "gemma_lmm"),
                run_test, args.genome_version,
                r["chrom"], r["pos"], r["variant_id"], r["original_id"],
                r["ref"], r["alt"], r["effect_allele"], r["other_allele"],
                fmt_num(r["beta"]), fmt_num(r["se"]), args.effect_scale,
                fmt_num(r["p_value"]), r["p_value_text"],
                fmt_num(r["neg_log10_p"]), fmt_num(r["eaf"]), fmt_num(r["maf"]),
                r["mac"] if r["mac"] is not None else "NA",
                r["n_obs"] if r["n_obs"] is not None else "NA",
                r["status"], r["filter_reason"],
                "true" if significant else "false",
            ])

    # ---- significant.tsv ----
    with open(out_dir / "significant.tsv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([c for c in GWAS_COLUMNS])
        for r in significant_rows:
            writer.writerow([
                SCHEMA_VERSION, args.analysis_id, args.trait,
                sample_meta.get("primary_model", "gemma_lmm"),
                run_test, args.genome_version,
                r["chrom"], r["pos"], r["variant_id"], r["original_id"],
                r["ref"], r["alt"], r["effect_allele"], r["other_allele"],
                fmt_num(r["beta"]), fmt_num(r["se"]), args.effect_scale,
                fmt_num(r["p_value"]), r["p_value_text"],
                fmt_num(r["neg_log10_p"]), fmt_num(r["eaf"]), fmt_num(r["maf"]),
                r["mac"] if r["mac"] is not None else "NA",
                r["n_obs"] if r["n_obs"] is not None else "NA",
                r["status"], r["filter_reason"], "true",
            ])

    # ---- invalid_variants.tsv ----
    with open(out_dir / "invalid_variants.tsv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["variant_id", "assoc_line", "status", "filter_reason"])
        for r in invalid_rows:
            writer.writerow([r["variant_id"], r["line_no"], r["status"], r["filter_reason"]])

    # ---- unreported_variants.tsv（待日志对账；不默认判定失败或通过）----
    with open(out_dir / "unreported_variants.tsv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["variant_id", "chrom", "pos", "reason_note"])
        for entry in sorted(unreported, key=lambda e: (e["chrom"], e["pos"])):
            writer.writerow([
                entry["variant_id"], entry["chrom"], entry["pos"],
                "未出现在 assoc 输出（待与软件日志/过滤记录对账；GEMMA -maf/-miss 过滤"
                "或性状级材料 MAC 不足均可能；不默认判定失败或通过）",
            ])

    # ---- metadata.json ----
    n_reported = len(harmonized)
    meta_out = {
        "schema": "gwasflow.harmonized.v1",
        "analysis_id": args.analysis_id,
        "trait": args.trait,
        "test": run_test,
        "p_column": args.p_column,
        "target_effect_allele": args.target_effect_allele,
        "status": status,
        "counts": {
            "n_reported": n_reported,
            "n_valid_p": m_tests,
            "n_invalid": n_invalid,
            "n_significant": len(significant_rows),
            "n_unreported": len(unreported),
            "n_total_variant_map": len(variant_map),
        },
        "multiple_testing": {
            "method": "Bonferroni（阈值 = alpha/(family_size×M)）",
            "alpha": str(alpha),
            "family_size": args.family_size,
            "m_tests": m_tests,
            "threshold": str(threshold) if threshold is not None else "NA",
            "note": "M=该性状有效 P 个数（记录实际检验数）；部分性状失败不得事后缩小 M",
        },
        "effect_scale": args.effect_scale,
        "genome_version": args.genome_version,
        "direction_conversion": {
            "policy": args.target_effect_allele,
            "n_converted": sum(1 for r in harmonized if r["converted"]),
            "rule": "beta×-1、eaf=1-eaf、等位基因对调、SE 与 P 不变",
        },
        "fingerprint": {
            "assoc_sha256": sha256_of_file(args.assoc),
            "variant_map_sha256": sha256_of_file(args.variant_map),
            "sample_metadata_sha256": sha256_of_file(args.sample_metadata),
            "run_metadata_sha256": sha256_of_file(args.run_metadata),
            "trait_mac_sha256": sha256_of_file(args.trait_mac),
        },
        "run_status_binding": run_status,
        "allow_partial_flag": bool(args.allow_partial),
        "schema_version_field": SCHEMA_VERSION,
    }
    with open(out_dir / "metadata.json", "w", encoding="utf-8") as handle:
        json.dump(meta_out, handle, indent=2, ensure_ascii=False)

    # ---- 退出码分级 ----
    print(
        f"[harmonize_gemma] 性状 '{args.trait}' 统一统计完成：reported={n_reported}，"
        f"有效 M={m_tests}，显著={len(significant_rows)}，无效={n_invalid}，"
        f"未报告={len(unreported)}，状态={status}"
    )
    if status == "partial":
        message = (
            f"存在 {n_invalid} 个数值异常行（invalid_variants.tsv；状态 partial）。"
            f"{'--allow-partial 生效，返回 0，但状态仍 partial。' if args.allow_partial else ''}"
        )
        print(f"[harmonize_gemma WARNING] {message}", file=sys.stderr)
        if not args.allow_partial:
            sys.exit(EXIT_PARTIAL)
    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
