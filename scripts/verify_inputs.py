#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GWASFlow 输入门控与文件校验（参数手册 §9.1 契约）。

产物：
  01.gate/input_gate.json      绑定生效配置摘要（SHA-256）与输入内容摘要（逐文件 MD5）
  01.gate/md5_check.tsv        MD5 三态（verified / FAIL / baseline_only）
  01.gate/reference_check.json 参考资源按需校验结果
  effective_config.yaml        Step 6 快照复制（由规则层移动到 13.summary/）

错误码：0=成功；2=格式/配置/输入错误（缺失/不一致/基线不符）。
红线：无基线输入绝不冒称 verified（baseline_only）。
"""

import argparse
import csv
import hashlib
import json
import os
import sys
from pathlib import Path

EXIT_OK = 0
EXIT_INPUT_ERROR = 2

TABLE_HEADERS = {
    "md5_check": ["path", "role", "expected_md5", "observed_md5", "status"],
}


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="verify_inputs.py",
        description=(
            "GWASFlow 输入门控：逐文件角色/组件/可读性校验、MD5 三态校验、"
            "表格契约校验与参考按需校验。"
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input-manifest", required=True,
                        help="输入清单 JSON（逐文件角色/路径/组件，由 snakefile 生成）")
    parser.add_argument("--samples", required=True, help="samples.csv 路径")
    parser.add_argument("--phenotypes", required=True, help="phenotypes.tsv 路径")
    parser.add_argument("--traits", required=True, help="traits.tsv 路径")
    parser.add_argument("--covariates", default=None, help="covariates.tsv 路径（性状声明协变量时必填）")
    parser.add_argument("--reference-manifest", required=True,
                        help="参考资源清单 JSON（组装/参考/注释摘要）")
    parser.add_argument("--expected-checksums", default=None,
                        help="提供方校验和基线 TSV（两列：path/md5 或 file/md5）")
    parser.add_argument("--checksum-mode", default="auto",
                        choices=["auto", "verified_required", "baseline"],
                        help="校验和模式：auto=有基线比对无基线标注；"
                             "verified_required=无基线阻断；baseline=仅计算不比对")
    parser.add_argument("--effective-config", required=True,
                        help="生效配置快照 YAML（Step 6 产物，绑定门控）")
    parser.add_argument("--threads", type=int, default=2, help="校验并行度")
    parser.add_argument("--required-covariates", default="false",
                        help="任一性状声明协变量时必须提供 covariates（true/false）")
    parser.add_argument("--out-dir", required=True, help="输出目录（必须为新目录或空目录）")
    return parser.parse_args(argv)


def fail(message, code=EXIT_INPUT_ERROR):
    print(f"[verify_inputs ERROR] {message}", file=sys.stderr)
    sys.exit(code)


def md5_of_file(path):
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_of_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def ensure_empty_dir(out_dir):
    path = Path(out_dir).resolve()
    if path.exists():
        if not path.is_dir():
            fail(f"输出路径已存在且不是目录: {path}")
        if any(path.iterdir()):
            fail(f"输出目录非空（拒绝覆盖既有产物）: {path}")
    else:
        path.mkdir(parents=True, exist_ok=True)
    return path


def load_json(path, label):
    p = Path(path)
    if not p.exists():
        fail(f"{label} 文件不存在: {p}")
    try:
        with open(p, encoding="utf-8") as handle:
            return json.load(handle)
    except json.JSONDecodeError as exc:
        fail(f"{label} JSON 解析失败: {exc}")


# ---- 逐文件角色与组件校验 ----

def check_file_entry(entry):
    """校验单个输入文件条目：角色/路径/组件/可读性。"""
    role = entry.get("role", "<未知角色>")
    path = entry.get("path")
    if not path:
        return [f"角色 {role} 缺少 path 字段"]
    p = Path(path)
    if not p.exists():
        return [f"角色 {role} 的文件不存在: {p}"]
    if not p.is_file():
        return [f"角色 {role} 的路径不是文件: {p}"]
    if not os.access(p, os.R_OK):
        return [f"角色 {role} 的文件不可读: {p}"]

    errors = []
    components = entry.get("components", [])
    fmt = entry.get("format", "")

    if fmt in ("vcf", "bcf") or role == "genotype_vcf":
        has_index = any((Path(str(path) + suf).exists()) for suf in (".csi", ".tbi"))
        if not has_index:
            errors.append(
                f"VCF/BCF 输入缺少索引（.csi 或 .tbi）: {path}"
            )
    if fmt == "bed" or role == "genotype_bed":
        prefix = entry.get("bfile") or str(path)
        for ext in (".bed", ".bim", ".fam"):
            if not Path(prefix + ext).exists():
                errors.append(f"BED 三件套不齐备，缺少 {prefix + ext}")
    return errors


# ---- MD5 三态 ----

def load_expected_checksums(path):
    """加载提供方校验和基线 TSV（两列无固定表头：第一列路径，第二列 MD5）。"""
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        fail(f"提供方校验和基线文件不存在: {p}")
    expected = {}
    with open(p, newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_no, row in enumerate(reader, start=1):
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 2:
                fail(f"校验和基线第 {line_no} 行列数不足: {row}")
            expected[row[0].strip()] = row[1].strip().lower()
    return expected


def run_md5_checks(files, expected, checksum_mode):
    """对每个输入文件执行 MD5 三态。

    verified: 基线存在且一致；FAIL: 基线存在但不一致（阻断）；
    baseline_only: 无基线（计算并标注，禁止冒称 verified）。
    """
    rows = []
    failures = []
    for entry in files:
        path = entry["path"]
        role = entry.get("role", "input")
        observed = md5_of_file(path)
        base = expected.get(path) or expected.get(os.path.basename(path))
        if checksum_mode == "baseline":
            status = "baseline_only"
        elif base is None:
            if checksum_mode == "verified_required":
                failures.append(f"verified_required 模式下缺少提供方基线: {path}")
                status = "FAIL"
            else:
                status = "baseline_only"
        elif base == observed:
            status = "verified"
        else:
            status = "FAIL"
            failures.append(
                f"MD5 校验失败: {path}（期望 {base}，实际 {observed}）"
            )
        rows.append({
            "path": path,
            "role": role,
            "expected_md5": base or "",
            "observed_md5": observed,
            "status": status,
        })
    return rows, failures


# ---- 表格契约校验 ----

def check_samples_table(path):
    errors = []
    seen = {}
    fid_iid = {}
    if not Path(path).exists():
        return [f"samples.csv 不存在: {path}"]
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = [c.strip() for c in (reader.fieldnames or [])]
        for col in ("sample", "fid", "iid"):
            if col not in fieldnames:
                errors.append(f"samples.csv 缺失必填列: {col}")
                return errors
        for line_no, row in enumerate(reader, start=2):
            row = {k: (v.strip() if v else "") for k, v in row.items()}
            sample = row.get("sample", "")
            if not sample:
                errors.append(f"samples.csv 第 {line_no} 行 sample 为空")
                continue
            if sample in seen:
                errors.append(
                    f"samples.csv sample 重复: '{sample}'（第 {line_no} 行与第 {seen[sample]} 行）"
                )
            seen[sample] = line_no
            pair = (row.get("fid", ""), row.get("iid", ""))
            if pair in fid_iid:
                errors.append(
                    f"samples.csv FID/IID 二元组重复: {pair}（第 {line_no} 行与第 {fid_iid[pair]} 行）"
                )
            fid_iid[pair] = line_no
    return errors


def check_unique_column(path, column, label):
    """校验 TSV 指定列唯一（错误含重复值与行号）。"""
    errors = []
    if not Path(path).exists():
        return [f"{label} 不存在: {path}"]
    seen = {}
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = [c.strip() for c in (reader.fieldnames or [])]
        if column not in fieldnames:
            return [f"{label} 缺失列: {column}（实际列: {fieldnames}）"]
        for line_no, row in enumerate(reader, start=2):
            value = (row.get(column) or "").strip()
            if not value:
                errors.append(f"{label} 第 {line_no} 行 {column} 为空")
                continue
            if value in seen:
                errors.append(
                    f"{label} {column} 重复: '{value}'（第 {line_no} 行与第 {seen[value]} 行）"
                )
            seen[value] = line_no
    return errors


# ---- 参考按需校验 ----

def check_references(reference_manifest, enabled_modules):
    """仅已启用模块所需参考资源参与校验（GFF3 在候选基因关闭时缺失不阻断）。"""
    version = reference_manifest.get("genome_version", "")
    resources = reference_manifest.get("resources", {})
    required_by_module = {
        "genotype_normalize": ["genome_fasta", "genome_fai", "chromosome_rename_map"],
        "candidate_genes": ["gff3"],
        "visualization": ["chromosome_table"],
    }
    required = set()
    for module in enabled_modules:
        required.update(required_by_module.get(module, []))
    # VCF/BCF 路线才需要 FASTA/FAI 与重命名映射；BED 路线 BIM 已是数值染色体
    if "genotype_normalize" in enabled_modules:
        fmt = reference_manifest.get("genotype_format", "vcf")
        if fmt not in ("vcf", "bcf"):
            required.discard("genome_fasta")
            required.discard("genome_fai")
            required.discard("chromosome_rename_map")

    report = {
        "genome_version": version,
        "genotype_format": reference_manifest.get("genotype_format", ""),
        "required_keys": sorted(required),
        "checked": {},
        "missing": [],
    }
    for key in sorted(required):
        path = resources.get(key)
        if not path:
            report["missing"].append({"key": key, "path": "<未配置>"})
            continue
        exists = Path(path).exists()
        entry = {"path": path, "exists": exists}
        if exists:
            entry["md5"] = md5_of_file(path)
            entry["size_bytes"] = os.path.getsize(path)
        report["checked"][key] = entry
        if not exists:
            report["missing"].append({"key": key, "path": path})
    return report


# ---- 主流程 ----

def main(argv=None):
    args = parse_args(argv)
    out_dir = ensure_empty_dir(args.out_dir)

    input_manifest = load_json(args.input_manifest, "输入清单")
    reference_manifest = load_json(args.reference_manifest, "参考资源清单")

    files = input_manifest.get("files", [])
    enabled_modules = input_manifest.get("enabled_modules", [])

    errors = []

    # 1. 逐文件角色/组件/可读性校验
    for entry in files:
        errors.extend(check_file_entry(entry))

    # 2. 表格契约校验（错误含重复值与行号）
    if args.required_covariates.lower() == "true" and not args.covariates:
        errors.append("任一性状声明协变量，必须提供 --covariates")
    table_paths = [
        ("samples.csv", args.samples),
        ("phenotypes.tsv", args.phenotypes),
        ("traits.tsv", args.traits),
    ]
    if args.covariates:
        table_paths.append(("covariates.tsv", args.covariates))

    errors.extend(check_samples_table(args.samples))
    errors.extend(check_unique_column(args.phenotypes, "sample", "phenotypes.tsv"))
    errors.extend(check_unique_column(args.traits, "trait", "traits.tsv"))
    if args.covariates:
        errors.extend(check_unique_column(args.covariates, "sample", "covariates.tsv"))

    if errors:
        for err in errors:
            print(f"[verify_inputs ERROR] {err}", file=sys.stderr)
        fail(f"输入校验失败（{len(errors)} 项错误）")

    # 3. MD5 三态
    expected = load_expected_checksums(args.expected_checksums)
    md5_rows, md5_failures = run_md5_checks(files, expected, args.checksum_mode)
    if md5_failures:
        for err in md5_failures:
            print(f"[verify_inputs ERROR] {err}", file=sys.stderr)
        fail("MD5 校验存在 FAIL 项")

    # 4. 参考按需校验
    reference_report = check_references(reference_manifest, enabled_modules)
    if reference_report["missing"]:
        for item in reference_report["missing"]:
            print(
                f"[verify_inputs ERROR] 参考资源缺失: {item['key']} -> {item['path']}",
                file=sys.stderr,
            )
        fail("参考资源按需校验失败")

    # 5. 门控产物绑定（生效配置摘要 + 输入内容摘要）
    config_sha256 = sha256_of_file(args.effective_config)
    gate = {
        "schema": "gwasflow.gate.v1",
        "config_sha256": config_sha256,
        "checksum_mode": args.checksum_mode,
        "inputs": [
            {
                "path": row["path"],
                "role": row["role"],
                "md5": row["observed_md5"],
                "size_bytes": os.path.getsize(row["path"]),
            }
            for row in md5_rows
        ],
        "tables": {
            label: {"path": str(path), "sha256": sha256_of_file(path)}
            for label, path in table_paths
        },
        "references": {
            "genome_version": reference_report["genome_version"],
            "checked_keys": sorted(reference_report["checked"].keys()),
        },
    }

    # 6. 产物落盘
    with open(out_dir / "input_gate.json", "w", encoding="utf-8") as handle:
        json.dump(gate, handle, indent=2, ensure_ascii=False)

    with open(out_dir / "md5_check.tsv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=TABLE_HEADERS["md5_check"], delimiter="\t"
        )
        writer.writeheader()
        for row in md5_rows:
            writer.writerow(row)

    with open(out_dir / "reference_check.json", "w", encoding="utf-8") as handle:
        json.dump(reference_report, handle, indent=2, ensure_ascii=False)

    # effective_config 快照复制（规则层移动到 13.summary/）
    with open(args.effective_config, encoding="utf-8") as src, open(
        out_dir / "effective_config.yaml", "w", encoding="utf-8"
    ) as dst:
        dst.write(src.read())

    n_verified = sum(1 for r in md5_rows if r["status"] == "verified")
    n_baseline = sum(1 for r in md5_rows if r["status"] == "baseline_only")
    print(
        f"[verify_inputs] 门控通过：{len(md5_rows)} 文件"
        f"（verified={n_verified}, baseline_only={n_baseline}），"
        f"产物目录 {out_dir}"
    )
    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
