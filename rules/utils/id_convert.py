#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GWASFlow 样本与性状定义加载（标准库实现，无 pandas 依赖）。

samples.csv 契约：sample,fid,iid 必填；可选 sample_name/group/include。
traits.tsv 契约（九列）：trait/phenotype_column/trait_type/unit/transform/
primary_model/numeric_covariates/categorical_covariates/n_pcs。
"""

import csv
import sys
from pathlib import Path

TRAIT_COLUMNS = [
    "trait",
    "phenotype_column",
    "trait_type",
    "unit",
    "transform",
    "primary_model",
    "numeric_covariates",
    "categorical_covariates",
    "n_pcs",
]

VALID_TRANSFORMS = {"none", "log", "log1p"}
VALID_TRAIT_TYPES = {"continuous"}
VALID_PRIMARY_MODELS = {"gemma_lmm"}


def _die(message: str, code: int = 2) -> None:
    print(f"[id_convert ERROR] {message}", file=sys.stderr)
    sys.exit(code)


def load_samples(csv_path, required_cols=None, index_col="sample"):
    """解析 samples.csv，返回以 sample 为键的字典。

    校验：必填列存在、sample 唯一、(fid,iid) 二元组唯一、include 布尔语义。
    错误信息含重复值与行号。
    """
    if required_cols is None:
        required_cols = ["sample", "fid", "iid"]

    file_path = Path(csv_path)
    if not file_path.exists():
        _die(f"找不到样本表文件: {file_path}")

    samples = {}
    fid_iid_seen = {}
    with open(file_path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        reader.fieldnames = [c.strip() for c in (reader.fieldnames or [])]

        missing = [c for c in required_cols if c not in reader.fieldnames]
        if missing:
            _die(f"样本表缺失必填列: {missing}（实际列: {reader.fieldnames}）")

        for line_no, row in enumerate(reader, start=2):
            row = {k: (v.strip() if v is not None else "") for k, v in row.items()}
            sample = row.get(index_col, "")
            if not sample:
                _die(f"样本表第 {line_no} 行 {index_col} 列为空")
            if sample in samples:
                _die(f"样本 ID 重复: '{sample}'（第 {line_no} 行，此前已出现）")

            fid = row.get("fid", "")
            iid = row.get("iid", "")
            if not fid or not iid:
                _die(f"样本 '{sample}'（第 {line_no} 行）fid/iid 不能为空")
            pair = (fid, iid)
            if pair in fid_iid_seen:
                prev = fid_iid_seen[pair]
                _die(
                    f"FID/IID 二元组重复: ('{fid}','{iid}') "
                    f"（第 {line_no} 行与第 {prev} 行）"
                )
            fid_iid_seen[pair] = line_no

            include_raw = row.get("include", "")
            if include_raw == "":
                include = True
            else:
                include = include_raw.lower() in {"true", "1", "yes"}
                if not include and include_raw.lower() not in {"false", "0", "no"}:
                    _die(
                        f"样本 '{sample}'（第 {line_no} 行）include 列值非法: "
                        f"'{include_raw}'（应为 true/false）"
                    )

            samples[sample] = {
                "sample": sample,
                "fid": fid,
                "iid": iid,
                "sample_name": row.get("sample_name", "") or sample,
                "group": row.get("group", ""),
                "include": include,
            }
    return samples


def load_traits(traits_tsv):
    """解析 traits.tsv（九列契约），返回以 trait ID 为键的定义字典。

    phenotype_column 缺省 = trait；transform 枚举校验；n_pcs 类型化；
    numeric/categorical_covariates 逗号列表解析。trait 重复即报错（含行号）。
    """
    file_path = Path(traits_tsv)
    if not file_path.exists():
        _die(f"找不到性状定义表: {file_path}")

    traits = {}
    with open(file_path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        reader.fieldnames = [c.strip() for c in (reader.fieldnames or [])]

        missing = [c for c in TRAIT_COLUMNS if c not in reader.fieldnames]
        if missing:
            _die(
                f"traits.tsv 缺失必填列: {missing}（实际列: {reader.fieldnames}）"
            )

        for line_no, row in enumerate(reader, start=2):
            row = {k: (v.strip() if v is not None else "") for k, v in row.items()}
            trait = row.get("trait", "")
            if not trait:
                _die(f"traits.tsv 第 {line_no} 行 trait 列为空")
            if trait in traits:
                _die(f"trait 重复: '{trait}'（第 {line_no} 行，此前已出现）")

            phenotype_column = row.get("phenotype_column") or trait
            transform = (row.get("transform") or "none").lower()
            if transform not in VALID_TRANSFORMS:
                _die(
                    f"trait '{trait}'（第 {line_no} 行）transform 非法: "
                    f"'{transform}'（合法值: {sorted(VALID_TRANSFORMS)}）"
                )
            trait_type = (row.get("trait_type") or "continuous").lower()
            if trait_type not in VALID_TRAIT_TYPES:
                _die(
                    f"trait '{trait}'（第 {line_no} 行）trait_type 非法: "
                    f"'{trait_type}'（首版仅支持 continuous）"
                )
            primary_model = (row.get("primary_model") or "gemma_lmm").lower()
            if primary_model not in VALID_PRIMARY_MODELS:
                _die(
                    f"trait '{trait}'（第 {line_no} 行）primary_model 非法: "
                    f"'{primary_model}'（首版仅支持 gemma_lmm）"
                )

            n_pcs_raw = row.get("n_pcs") or "0"
            try:
                n_pcs = int(n_pcs_raw)
            except ValueError:
                _die(
                    f"trait '{trait}'（第 {line_no} 行）n_pcs 非整数: '{n_pcs_raw}'"
                )
            if n_pcs < 0:
                _die(f"trait '{trait}'（第 {line_no} 行）n_pcs 不能为负: {n_pcs}")

            def _split_list(field):
                raw = row.get(field) or ""
                return [x.strip() for x in raw.split(",") if x.strip()]

            numeric_covariates = _split_list("numeric_covariates")
            categorical_covariates = _split_list("categorical_covariates")
            overlap = set(numeric_covariates) & set(categorical_covariates)
            if overlap:
                _die(
                    f"trait '{trait}'（第 {line_no} 行）数值与分类协变量重复: "
                    f"{sorted(overlap)}"
                )

            unit = row.get("unit") or "unspecified"
            if transform != "none":
                effect_scale = f"{transform}({unit})"
            else:
                effect_scale = unit

            traits[trait] = {
                "trait": trait,
                "phenotype_column": phenotype_column,
                "trait_type": trait_type,
                "unit": unit,
                "transform": transform,
                "primary_model": primary_model,
                "numeric_covariates": numeric_covariates,
                "categorical_covariates": categorical_covariates,
                "n_pcs": n_pcs,
                "effect_scale": effect_scale,
            }
    return traits
