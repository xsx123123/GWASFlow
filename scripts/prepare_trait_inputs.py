#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GWASFlow 性状级输入准备（参数手册 §2 契约重实现，M1）。

以最终 QC 的 FAM 顺序为唯一对齐权威，全矩阵+缺失掩码：
- 表型/协变量/PC 全部按 FAM 顺序输出完整 N 行（禁止截短表型却传原 N 的 BED/K）
- 无效材料保留行位置、表型置 NA；有效材料生成 keep.tsv
- kinship_all.pheno.txt 全 1（全队列 K）

红线：
- -9 默认合法负值（显式加入 --missing-values 才视为缺失）
- 变换值域校验（log: y>0；log1p: y>-1，违反退出 2）
- 分类协变量 K-1 dummy（纳入材料内字典序参考组）
- 设计矩阵共线/秩不足/完整 N<=列数+1/残余自由度不足 → 退出 2 并提示共线来源列
- PC 身份精确匹配（含 FID 二元组；无 FID 时要求 FAM 中 IID 唯一）
- 有效材料 < min-samples → 退出 2（含当前 N 与阈值与缺失占比）

错误码：0=成功；2=格式/配置/输入错误。
"""

import argparse
import csv
import hashlib
import json
import math
import os
import sys
from pathlib import Path

EXIT_OK = 0
EXIT_INPUT_ERROR = 2
DEFAULT_MISSING = ["NA", "NAN", "."]


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="prepare_trait_inputs.py",
        description=(
            "以最终 QC 的 FAM 顺序为唯一权威，对齐材料级数量性状、协变量与 PCs，"
            "建立缺失掩码与 GEMMA 输入。不读 BED，不计算 PCA，不计算 kinship。"
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--fam", required=True, help="最终 QC .fam（六列；FID/IID 二元组唯一）")
    parser.add_argument("--samples", required=True, help="samples.csv（sample,fid,iid）")
    parser.add_argument("--phenotypes", required=True, help="phenotypes.tsv（sample 和表型列）")
    parser.add_argument("--trait", required=True, help="性状内部 ID（写入 metadata）")
    parser.add_argument("--phenotype-column", default=None, help="实际表型列名（默认同 trait）")
    parser.add_argument("--covariates", default=None, help="covariates.tsv")
    parser.add_argument("--numeric-covariates", default="", help="连续数值列（逗号列表）")
    parser.add_argument("--categorical-covariates", default="", help="分类列（逗号列表）")
    parser.add_argument("--pcs", default=None, help="PCA eigenvec 文件（接受 #FID/IID 或 #IID 表头）")
    parser.add_argument("--n-pcs", type=int, default=0, help="选择 PC1..PCn（不能为负）")
    parser.add_argument("--skip-pcs", action="store_true",
                        help="only_qc 降级：跳过 PC（metadata 记录 pc_skipped）")
    parser.add_argument("--transform", default="none", choices=["none", "log", "log1p"],
                        help="表型变换（log 要求 y>0；log1p 要求 y>-1）")
    parser.add_argument("--unit", default="unspecified", help="原始性状单位")
    parser.add_argument("--missing-values", default="NA,NaN,.",
                        help="缺失值标记（逗号列表，大小写不敏感；空值始终视为缺失）")
    parser.add_argument("--min-samples", type=int, default=20,
                        help="最少完整材料（>=3；文档注明不代表足够统计功效）")
    parser.add_argument("--rank-tolerance", type=float, default=1e-10,
                        help="归一化列的数值秩判定容差（0<值<1）")
    parser.add_argument("--out-dir", required=True, help="输出根目录（必须为新目录或空目录）")
    return parser.parse_args(argv)


def fail(message, code=EXIT_INPUT_ERROR):
    print(f"[prepare_trait_inputs ERROR] {message}", file=sys.stderr)
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


def is_missing(value, missing_set):
    if value is None:
        return True
    text = str(value).strip()
    if text == "":
        return True
    return text.upper() in missing_set


def parse_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def read_fam(fam_path):
    """读 FAM（顺序权威）：[(fid, iid), ...]，同时校验二元组唯一。"""
    entries = []
    seen = {}
    with open(fam_path, encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            parts = line.split()
            if len(parts) < 2:
                fail(f"FAM 第 {line_no} 行列数不足: {line.strip()!r}")
            pair = (parts[0], parts[1])
            if pair in seen:
                fail(f"FAM FID/IID 二元组重复: {pair}（第 {line_no} 行与第 {seen[pair]} 行）")
            seen[pair] = line_no
            entries.append(pair)
    return entries


def read_samples(csv_path):
    """读 samples.csv：{sample: {fid, iid, include, ...}} 与 (fid,iid)→sample 反查。"""
    samples = {}
    fid_iid_map = {}
    with open(csv_path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for line_no, row in enumerate(reader, start=2):
            row = {k: (v.strip() if v else "") for k, v in row.items()}
            sample = row.get("sample", "")
            if not sample:
                fail(f"samples.csv 第 {line_no} 行 sample 为空")
            fid, iid = row.get("fid", ""), row.get("iid", "")
            pair = (fid, iid)
            if sample in samples:
                fail(f"samples.csv sample 重复: '{sample}'（第 {line_no} 行）")
            if pair in fid_iid_map:
                fail(
                    f"samples.csv FID/IID 二元组重复: {pair}"
                    f"（第 {line_no} 行与第 {fid_iid_map[pair]} 行）"
                )
            include_raw = row.get("include", "")
            include = True if include_raw == "" else include_raw.lower() in {"true", "1", "yes"}
            samples[sample] = {"fid": fid, "iid": iid, "include": include}
            fid_iid_map[pair] = sample
    return samples, fid_iid_map


def read_phenotypes(tsv_path, sample_registry):
    """读 phenotypes.tsv：{sample: {column: raw_value}}；未知 sample 报错。"""
    table = {}
    with open(tsv_path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = [c.strip() for c in (reader.fieldnames or [])]
        if "sample" not in fieldnames:
            fail(f"phenotypes.tsv 缺失 sample 列（实际列: {fieldnames}）")
        for line_no, row in enumerate(reader, start=2):
            row = {k: (v.strip() if v else "") for k, v in row.items()}
            sample = row.get("sample", "")
            if sample in table:
                fail(f"phenotypes.tsv sample 重复: '{sample}'（第 {line_no} 行）")
            if sample not in sample_registry:
                fail(
                    f"phenotypes.tsv 第 {line_no} 行出现未知 sample: '{sample}'"
                    f"（不在 samples.csv 登记中）"
                )
            table[sample] = row
    return table, fieldnames


def read_covariates(tsv_path, sample_registry):
    table = {}
    with open(tsv_path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = [c.strip() for c in (reader.fieldnames or [])]
        if "sample" not in fieldnames:
            fail(f"covariates.tsv 缺失 sample 列（实际列: {fieldnames}）")
        for line_no, row in enumerate(reader, start=2):
            row = {k: (v.strip() if v else "") for k, v in row.items()}
            sample = row.get("sample", "")
            if sample in table:
                fail(f"covariates.tsv sample 重复: '{sample}'（第 {line_no} 行）")
            if sample not in sample_registry:
                fail(
                    f"covariates.tsv 第 {line_no} 行出现未知 sample: '{sample}'"
                    f"（不在 samples.csv 登记中；QC 外材料允许出现，但必须登记）"
                )
            table[sample] = row
    return table, fieldnames


def read_pcs(pcs_path, n_pcs):
    """读 PCA eigenvec：{（fid,iid）或 iid: [pc1..pcn]}。

    含 FID 表头（#FID IID PC1...）→ 二元组键；无 FID（#IID PC1...）→ iid 键。
    """
    result = {}
    has_fid = False
    with open(pcs_path, encoding="utf-8") as handle:
        header = handle.readline()
        cols = [c for c in header.strip().lstrip("#").split() if c]
        if len(cols) >= 2 and cols[0].upper() == "FID" and cols[1].upper() == "IID":
            has_fid = True
            pc_start = 2
        elif cols and cols[0].upper() == "IID":
            has_fid = False
            pc_start = 1
        else:
            fail(f"PCA eigenvec 表头无法解析（期望 #FID IID PC... 或 #IID PC...）: {header!r}")
        available_pcs = len(cols) - pc_start
        if available_pcs < n_pcs:
            fail(
                f"PCA eigenvec 仅含 {available_pcs} 个 PC，少于请求的 n_pcs={n_pcs}"
            )
        for line_no, line in enumerate(handle, start=2):
            parts = line.split()
            if not parts:
                continue
            if has_fid:
                if len(parts) < pc_start + n_pcs:
                    fail(f"eigenvec 第 {line_no} 行 PC 列不足")
                key = (parts[0], parts[1])
                values = [parse_float(v) for v in parts[pc_start:pc_start + n_pcs]]
            else:
                if len(parts) < pc_start + n_pcs:
                    fail(f"eigenvec 第 {line_no} 行 PC 列不足")
                key = parts[0]
                values = [parse_float(v) for v in parts[pc_start:pc_start + n_pcs]]
            if key in result:
                fail(f"eigenvec 第 {line_no} 行身份重复: {key}")
            if any(v is None for v in values):
                fail(f"eigenvec 第 {line_no} 行 PC 值非数值")
            result[key] = values
    return result, has_fid


# ---- 秩检查（标准库实现，归一化列 + 高斯消元数值秩）----

def matrix_rank(columns, tolerance):
    """列向量矩阵的数值秩（归一化后高斯消元）。

    columns: list[list[float]]，等长。
    返回 (rank, collinear_groups)——共线组用于错误信息定位来源列。
    """
    if not columns:
        return 0, []
    n_rows = len(columns[0])
    # 归一化每列
    normalized = []
    for col in columns:
        norm = math.sqrt(sum(x * x for x in col))
        if norm <= tolerance * max(1.0, abs(max(col, key=abs, default=0.0))):
            normalized.append([0.0] * n_rows)
        else:
            normalized.append([x / norm for x in col])
    # 高斯消元（对转置矩阵：行=原始列）
    basis = []
    collinear = []
    for idx, vec in enumerate(normalized):
        v = list(vec)
        for b_idx, b in enumerate(basis):
            dot = sum(x * y for x, y in zip(v, b))
            v = [x - dot * y for x, y in zip(v, b)]
        norm = math.sqrt(sum(x * x for x in v))
        if norm > tolerance:
            basis.append([x / norm for x in v])
        else:
            collinear.append(idx)
    return len(basis), collinear


def main(argv=None):
    args = parse_args(argv)
    out_dir = ensure_empty_dir(args.out_dir)

    if args.min_samples < 3:
        fail(f"--min-samples 必须 >=3（当前 {args.min_samples}；该阈值不代表足够统计功效）")
    if not (0 < args.rank_tolerance < 1):
        fail(f"--rank-tolerance 必须在 (0,1) 区间（当前 {args.rank_tolerance}）")
    if args.n_pcs < 0:
        fail(f"--n-pcs 不能为负（当前 {args.n_pcs}）")

    phenotype_column = args.phenotype_column or args.trait
    missing_set = {m.strip().upper() for m in args.missing_values.split(",") if m.strip()}
    # 空值始终视为缺失（is_missing 已处理）；-9 不在默认缺失集（合法负值）

    # ---- 数据加载 ----
    fam_order = read_fam(args.fam)
    samples_registry, fid_iid_to_sample = read_samples(args.samples)
    phenotypes, pheno_columns = read_phenotypes(args.phenotypes, samples_registry)

    if phenotype_column not in pheno_columns:
        fail(
            f"phenotypes.tsv 缺失表型列 '{phenotype_column}'"
            f"（可用列: {pheno_columns}）"
        )

    covariates, cov_columns = ({}, [])
    if args.covariates:
        covariates, cov_columns = read_covariates(args.covariates, samples_registry)

    numeric_cov_names = [c.strip() for c in args.numeric_covariates.split(",") if c.strip()]
    categorical_cov_names = [
        c.strip() for c in args.categorical_covariates.split(",") if c.strip()
    ]
    overlap = set(numeric_cov_names) & set(categorical_cov_names)
    if overlap:
        fail(f"数值与分类协变量重复: {sorted(overlap)}")

    pc_names = [f"PC{i + 1}" for i in range(args.n_pcs)]
    name_conflicts = set(numeric_cov_names) | set(categorical_cov_names)
    conflicts = name_conflicts & set(pc_names)
    if conflicts:
        fail(f"协变量列与 PC 名重复: {sorted(conflicts)}")

    for col in numeric_cov_names + categorical_cov_names:
        if col not in cov_columns:
            fail(f"covariates.tsv 缺失声明列 '{col}'（可用列: {cov_columns}）")
        if col == "sample":
            fail("协变量列不能使用保留列名 'sample'")

    pcs_data, pcs_has_fid = ({}, False)
    pc_skipped = False
    if args.skip_pcs:
        pc_skipped = True
        if args.n_pcs > 0:
            print(
                "[prepare_trait_inputs WARNING] --skip-pcs 生效："
                "跳过 PC 协变量（only_qc 降级路径，metadata 记录 pc_skipped）",
                file=sys.stderr,
            )
    elif args.n_pcs > 0:
        if not args.pcs:
            fail(f"n_pcs={args.n_pcs} 但未提供 --pcs（或使用 --skip-pcs 降级）")
        pcs_data, pcs_has_fid = read_pcs(args.pcs, args.n_pcs)
        if not pcs_has_fid:
            iid_counts = {}
            for fid, iid in fam_order:
                iid_counts[iid] = iid_counts.get(iid, 0) + 1
            duplicated = [k for k, v in iid_counts.items() if v > 1]
            if duplicated:
                fail(
                    f"PCA eigenvec 无 FID 且 FAM 中 IID 不唯一: {duplicated}"
                    f"（无法安全匹配 PC）"
                )

    # ---- 逐 FAM 材料构建（全矩阵 + 缺失掩码）----
    rows = []
    for fid, iid in fam_order:
        sample = fid_iid_to_sample.get((fid, iid))
        entry = {
            "fid": fid, "iid": iid, "sample": sample or "",
            "included": True, "exclusion_reason": "",
            "phenotype_raw": None, "phenotype_value": None,
            "numeric_covs": {}, "categorical_covs": {}, "pcs": [],
        }
        if sample is None:
            entry["included"] = False
            entry["exclusion_reason"] = "fam_identity_not_in_samples_table"
        elif not samples_registry[sample]["include"]:
            entry["included"] = False
            entry["exclusion_reason"] = "excluded_by_sample_sheet(include=false)"

        if entry["included"]:
            pheno_row = phenotypes.get(sample)
            if pheno_row is None:
                entry["included"] = False
                entry["exclusion_reason"] = "phenotype_missing(sample_not_in_table)"
            else:
                raw = pheno_row.get(phenotype_column, "")
                if is_missing(raw, missing_set):
                    entry["included"] = False
                    entry["exclusion_reason"] = "phenotype_missing"
                else:
                    value = parse_float(raw)
                    if value is None:
                        entry["included"] = False
                        entry["exclusion_reason"] = f"phenotype_non_numeric({raw!r})"
                    else:
                        entry["phenotype_raw"] = raw
                        entry["phenotype_value"] = value

        if entry["included"] and args.covariates:
            cov_row = covariates.get(sample)
            if cov_row is None:
                entry["included"] = False
                entry["exclusion_reason"] = "covariate_missing(sample_not_in_table)"
            else:
                for col in numeric_cov_names:
                    raw = cov_row.get(col, "")
                    if is_missing(raw, missing_set):
                        entry["included"] = False
                        entry["exclusion_reason"] = f"covariate_missing:{col}"
                        break
                    value = parse_float(raw)
                    if value is None:
                        entry["included"] = False
                        entry["exclusion_reason"] = f"covariate_non_numeric:{col}({raw!r})"
                        break
                    entry["numeric_covs"][col] = value

        if entry["included"] and categorical_cov_names:
            cov_row = covariates.get(sample) if args.covariates else None
            if cov_row is None and categorical_cov_names:
                entry["included"] = False
                entry["exclusion_reason"] = "covariate_missing(covariates_not_provided)"
            else:
                for col in categorical_cov_names:
                    raw = (cov_row or {}).get(col, "")
                    if is_missing(raw, missing_set):
                        entry["included"] = False
                        entry["exclusion_reason"] = f"covariate_missing:{col}"
                        break
                    entry["categorical_covs"][col] = raw

        if entry["included"] and args.n_pcs > 0 and not pc_skipped:
            key = (fid, iid) if pcs_has_fid else iid
            pc_values = pcs_data.get(key)
            if pc_values is None:
                entry["included"] = False
                entry["exclusion_reason"] = "pc_missing(identity_not_in_eigenvec)"
            else:
                entry["pcs"] = pc_values

        rows.append(entry)

    # ---- 变换值域校验 ----
    valid_rows = [r for r in rows if r["included"]]
    if args.transform == "log":
        bad = [r for r in valid_rows if r["phenotype_value"] <= 0]
        if bad:
            fail(
                f"transform=log 要求表型值 > 0，违反材料数 {len(bad)}"
                f"（首个: {bad[0]['fid']}/{bad[0]['iid']} 值 {bad[0]['phenotype_value']}）"
            )
    elif args.transform == "log1p":
        bad = [r for r in valid_rows if r["phenotype_value"] <= -1]
        if bad:
            fail(
                f"transform=log1p 要求表型值 > -1，违反材料数 {len(bad)}"
                f"（首个: {bad[0]['fid']}/{bad[0]['iid']} 值 {bad[0]['phenotype_value']}）"
            )

    # ---- 表型零方差阻断 ----
    if valid_rows:
        values = [r["phenotype_value"] for r in valid_rows]
        if len(set(values)) == 1:
            fail(
                f"表型零方差（全部有效材料值 = {values[0]}），模型不可估计"
            )
    else:
        fail("有效材料为 0（全部被排除，详见 sample_inclusion 排除原因）")

    # ---- 分类协变量编码（纳入材料内字典序参考组，K-1 dummy）----
    n_valid = len(valid_rows)
    design_columns = ["intercept"]
    encoding = {}
    dummy_columns = {}
    for col in categorical_cov_names:
        levels = sorted({r["categorical_covs"][col] for r in valid_rows})
        if len(levels) < 2:
            fail(
                f"分类协变量 '{col}' 在纳入材料内仅 {len(levels)} 个有效水平（{levels}），"
                f"无法 K-1 编码"
            )
        reference = levels[0]
        encoding[col] = {"reference_level": reference, "levels": levels,
                         "method": "K-1 dummy（字典序参考组；需固定参考组请见文档声明字段）"}
        for level in levels[1:]:
            col_name = f"{col}[{level}]"
            dummy_columns[col_name] = (col, level)
            design_columns.append(col_name)
    for col in numeric_cov_names:
        design_columns.append(col)
    design_columns.extend(pc_names)

    # ---- 设计合法性检查 ----
    n_design_cols = len(design_columns) - 1   # 除截距
    if n_valid <= n_design_cols + 1:
        fail(
            f"完整有效 N（{n_valid}）<= 协变量列数+1（{n_design_cols + 1}），"
            f"残余自由度不足（设计列: {design_columns}）"
        )

    # 秩检查（数值协变量 + dummy + PC；共线来源列定位）
    design_matrix_cols = []
    design_matrix_names = []
    for name in design_columns[1:]:
        if name in dummy_columns:
            col, level = dummy_columns[name]
            design_matrix_cols.append(
                [1.0 if r["categorical_covs"][col] == level else 0.0 for r in valid_rows]
            )
        elif name in numeric_cov_names:
            design_matrix_cols.append([r["numeric_covs"][name] for r in valid_rows])
        else:
            pc_idx = int(name[2:]) - 1
            design_matrix_cols.append([r["pcs"][pc_idx] for r in valid_rows])
        design_matrix_names.append(name)

    if design_matrix_cols:
        rank, collinear_idx = matrix_rank(design_matrix_cols, args.rank_tolerance)
        if rank < len(design_matrix_cols):
            sources = [design_matrix_names[i] for i in collinear_idx]
            fail(
                f"设计矩阵共线（数值秩 {rank} < 列数 {len(design_matrix_cols)}），"
                f"共线来源列: {sources}（含截距外设计列 {design_matrix_names}）"
            )

    # ---- min-samples 门槛 ----
    if n_valid < args.min_samples:
        n_total = len(rows)
        fail(
            f"有效材料 {n_valid} < min_samples {args.min_samples}"
            f"（FAM 总数 {n_total}，缺失/排除占比 "
            f"{(n_total - n_valid) / n_total:.1%}；注意该阈值不代表足够统计功效）"
        )

    # ---- 产物写入（八件）----
    def fmt(value):
        return repr(value) if value == value else "NA"

    # 1. phenotype.gemma.txt（完整 N 行；无效 NA；应用变换）
    with open(out_dir / "phenotype.gemma.txt", "w", encoding="utf-8") as handle:
        for r in rows:
            if r["included"]:
                v = r["phenotype_value"]
                if args.transform == "log":
                    v = math.log(v)
                elif args.transform == "log1p":
                    v = math.log1p(v)
                handle.write(f"{v!r}\n")
            else:
                handle.write("NA\n")

    # 2. covariates.gemma.txt（完整 N 行；第一列截距；无效材料全 NA 行为：仅截距+NA）
    with open(out_dir / "covariates.gemma.txt", "w", encoding="utf-8") as handle:
        for r in rows:
            parts = ["1"]
            if r["included"]:
                for name in design_columns[1:]:
                    if name in dummy_columns:
                        col, level = dummy_columns[name]
                        parts.append("1" if r["categorical_covs"][col] == level else "0")
                    elif name in numeric_cov_names:
                        parts.append(f"{r['numeric_covs'][name]!r}")
                    else:
                        pc_idx = int(name[2:]) - 1
                        parts.append(f"{r['pcs'][pc_idx]!r}")
            else:
                # 无效材料：协变量置 NA（保留行位置；GEMMA 按表型 NA 掩码）
                parts.extend(["NA"] * (len(design_columns) - 1))
            handle.write("\t".join(parts) + "\n")

    # 3. kinship_all.pheno.txt（完整 N 行 1，全队列 K）
    with open(out_dir / "kinship_all.pheno.txt", "w", encoding="utf-8") as handle:
        for _ in rows:
            handle.write("1\n")

    # 4. keep.tsv（有效 FID/IID，无表头）
    with open(out_dir / "keep.tsv", "w", encoding="utf-8") as handle:
        for r in valid_rows:
            handle.write(f"{r['fid']}\t{r['iid']}\n")

    # 5. sample_inclusion.tsv（全部 FAM 材料 + 逐项排除原因；行序与 FAM 一致）
    with open(out_dir / "sample_inclusion.tsv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["fid", "iid", "sample", "included", "exclusion_reason"])
        for r in rows:
            writer.writerow([
                r["fid"], r["iid"], r["sample"],
                "true" if r["included"] else "false",
                r["exclusion_reason"],
            ])

    # 6. sample_order.tsv（完整 FAM 顺序，有表头）
    order_text = "".join(f"{r['fid']}\t{r['iid']}\n" for r in rows)
    with open(out_dir / "sample_order.tsv", "w", encoding="utf-8") as handle:
        handle.write("fid\tiid\n")
        handle.write(order_text)

    # 7. design_columns.tsv（协变量矩阵列名和顺序）
    with open(out_dir / "design_columns.tsv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["column", "type", "source"])
        writer.writerow(["intercept", "intercept", "explicit"])
        for name in design_columns[1:]:
            if name in dummy_columns:
                col, _ = dummy_columns[name]
                writer.writerow([name, "categorical_dummy", f"covariates.tsv:{col}"])
            elif name in numeric_cov_names:
                writer.writerow([name, "numeric", f"covariates.tsv:{name}"])
            else:
                writer.writerow([name, "pc", "pca.eigenvec"])

    # 8. metadata.json
    effect_scale = (
        f"{args.transform}({args.unit})" if args.transform != "none" else args.unit
    )
    order_sha = hashlib.sha256(order_text.encode()).hexdigest()
    sources = {
        "fam": {"path": str(Path(args.fam).resolve()), "sha256": sha256_of_file(args.fam)},
        "samples": {"path": str(Path(args.samples).resolve()),
                    "sha256": sha256_of_file(args.samples)},
        "phenotypes": {"path": str(Path(args.phenotypes).resolve()),
                       "sha256": sha256_of_file(args.phenotypes)},
    }
    if args.covariates:
        sources["covariates"] = {
            "path": str(Path(args.covariates).resolve()),
            "sha256": sha256_of_file(args.covariates),
        }
    if args.pcs and args.n_pcs > 0 and not pc_skipped:
        sources["pcs"] = {"path": str(Path(args.pcs).resolve()),
                          "sha256": sha256_of_file(args.pcs)}

    metadata = {
        "schema": "gwasflow.trait_prep.v1",
        "trait": args.trait,
        "phenotype_column": phenotype_column,
        "transform": args.transform,
        "unit": args.unit,
        "effect_scale": effect_scale,
        "n_total": len(rows),
        "n_valid": n_valid,
        "sample_order_sha256": order_sha,
        "sample_order_n": len(rows),
        "categorical_encoding": encoding,
        "design_columns": design_columns,
        "n_pcs": 0 if pc_skipped else args.n_pcs,
        "pc_skipped": pc_skipped,
        "missing_values": sorted(missing_set),
        "min_samples": args.min_samples,
        "rank_tolerance": args.rank_tolerance,
        "fam_authority": "FAM 顺序为唯一对齐权威（全矩阵+缺失掩码路线）",
        "sources": sources,
    }
    with open(out_dir / "metadata.json", "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, ensure_ascii=False)

    # 产物 SHA-256 追记
    products = {}
    for name in (
        "phenotype.gemma.txt", "covariates.gemma.txt", "kinship_all.pheno.txt",
        "keep.tsv", "sample_inclusion.tsv", "sample_order.tsv",
        "design_columns.tsv",
    ):
        products[name] = sha256_of_file(out_dir / name)
    metadata["products"] = products
    with open(out_dir / "metadata.json", "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, ensure_ascii=False)

    print(
        f"[prepare_trait_inputs] 性状 '{args.trait}' 准备完成："
        f"有效 {n_valid}/{len(rows)} 材料，产物目录 {out_dir}"
    )
    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
