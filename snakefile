#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# GWASFlow 入口（小写 snakefile）：九步初始化
# 框架对齐：FlowFrame v2.0（CygnusX-Labs/FlowFrame, commit 24a0cd3e）
# 版本唯一来源：config/config.yaml 的 pipeline_version（本文件不写版本号字符串）

import os
import sys

from snakemake.utils import min_version, validate

# ---- Step 1: 导入框架工具 ----
from rules.utils.id_convert import load_samples, load_traits
from rules.utils.validate import (
    load_user_config,
    validate_genome_version,
    validate_module_switches,
)
from rules.utils.reference_update import (
    resolve_reference_paths,
    check_gwas_reference_paths,
)
from rules.utils.resource_manager import rule_resource
from rules.utils.gwas_targets import (
    build_target_plan,
    resolve_modules,
    format_module_status,
)
from rules.utils.common import DataDeliver

# Lock Snakemake Version
min_version("9.9.0")

# ---- Logger 三级降级链（D9：禁止顶层裸 import 插件）----
try:
    from snakemake_logger_plugin_rich_loguru import get_analysis_logger
    logger = get_analysis_logger()
except ImportError:
    try:
        from loguru import logger
    except ImportError:
        import logging

        logging.basicConfig(level="INFO")
        logger = logging.getLogger("GWASFlow")

# ---- Step 2: 加载默认四件套（声明顺序 = 优先级低 → 高）----
configfile: "config/config.yaml"
configfile: "config/reference.yaml"
configfile: "config/run_parameter.yaml"
configfile: "config/cluster_config.yaml"

# ---- Step 3: analysisyaml 覆盖（最高优先级）----
load_user_config(config, cmd_arg_name="analysisyaml")

# ---- Step 4: 解析参考路径（GWAS_reference 相对路径 → 绝对路径）----
resolve_reference_paths(config, index_keys=["GWAS_reference"])

# ---- Step 5: Schema 验证 + Genome_Version 校验 + 未实现模块开关拒绝 ----
validate(config, "schema/config.schema.yaml")
validate_genome_version(config=config, logger=logger)
validate_module_switches(config=config, logger=logger)

# ---- Step 6: 生效配置快照与模块开关解析输出 ----
# 深合并终值序列化（落盘延迟到工作区 02 门控规则 output，避免仓库目录副作用）
import yaml as _yaml

EFFECTIVE_CONFIG_TEXT = _yaml.safe_dump(
    config, default_flow_style=False, allow_unicode=True, sort_keys=True
)

MODULE_STATUS_PREVIEW = resolve_modules(config, use_pcs=None)

# 生效开关清单与自动启用原因（控制台/日志，spec §4.4-3/4）
logger.info("GWASFlow 生效模块开关清单：\n" + format_module_status(MODULE_STATUS_PREVIEW))

# 快照源文件写入 workflow 绝对路径（workdir 重定向前；02 规则以 input 挂接）
import json as _json

EFFECTIVE_CONFIG_SOURCE = os.path.join(
    config["workflow"], "13.summary", ".effective_config.source.yaml"
)
os.makedirs(os.path.dirname(EFFECTIVE_CONFIG_SOURCE), exist_ok=True)
with open(EFFECTIVE_CONFIG_SOURCE, "w", encoding="utf-8") as handle:
    handle.write(EFFECTIVE_CONFIG_TEXT)


def _build_input_manifest(cfg, modules):
    fmt = cfg.get("genotype_format", "vcf")
    files = []
    if fmt in ("vcf", "bcf"):
        files.append({
            "role": "genotype",
            "format": fmt,
            "path": cfg.get("genotype_input", ""),
            "components": ["vcf", "index"],
        })
    else:
        files.append({
            "role": "genotype",
            "format": "bed",
            "path": (cfg.get("genotype_bfile") or "") + ".bed",
            "bfile": cfg.get("genotype_bfile", ""),
            "components": ["bed", "bim", "fam"],
        })
        if cfg.get("input_variant_map"):
            files.append({
                "role": "variant_map",
                "format": "tsv",
                "path": cfg["input_variant_map"],
                "components": ["tsv"],
            })
    for role, key in (
        ("samples", "sample_csv"),
        ("phenotypes", "phenotypes_tsv"),
        ("traits", "traits_tsv"),
        ("covariates", "covariates_tsv"),
    ):
        if cfg.get(key):
            files.append({"role": role, "format": "table", "path": cfg[key]})
    if cfg.get("expected_checksums_tsv"):
        files.append({"role": "expected_checksums", "format": "tsv",
                      "path": cfg["expected_checksums_tsv"]})
    return {
        "schema": "gwasflow.input_manifest.v1",
        "enabled_modules": sorted(m for m, s in modules.items() if s.enabled),
        "files": files,
    }


INPUT_MANIFEST_SOURCE = os.path.join(
    config["workflow"], "13.summary", ".input_manifest.source.json"
)
INPUT_MANIFEST = _build_input_manifest(config, MODULE_STATUS_PREVIEW)
with open(INPUT_MANIFEST_SOURCE, "w", encoding="utf-8") as handle:
    _json.dump(INPUT_MANIFEST, handle, indent=2, ensure_ascii=False)

# 参考资源清单（GWAS_reference 当前版本 + 基因型格式）
REFERENCE_MANIFEST_SOURCE = os.path.join(
    config["workflow"], "13.summary", ".reference_manifest.source.json"
)
_REFERENCE_VERSION = str(config.get("Genome_Version", "")).strip()
_REFERENCE_MANIFEST = {
    "schema": "gwasflow.reference_manifest.v1",
    "genome_version": _REFERENCE_VERSION,
    "genotype_format": config.get("genotype_format", "vcf"),
    "resources": (config.get("GWAS_reference", {}) or {}).get(_REFERENCE_VERSION, {}),
}
with open(REFERENCE_MANIFEST_SOURCE, "w", encoding="utf-8") as handle:
    _json.dump(_REFERENCE_MANIFEST, handle, indent=2, ensure_ascii=False)

# ---- Step 7: 工作目录重定向（不可逆全局副作用；此后相对路径按新工作目录解析）----
workdir: config["workflow"]
logger.info(f"Redirect workspaces to {config['workflow']}")

# ---- Step 8: 加载样本与性状（解析发生在 workdir 重定向后，输入表为绝对路径）----
SAMPLES = load_samples(
    config["sample_csv"], required_cols=["sample", "fid", "iid"]
)
TRAITS = load_traits(config["traits_tsv"])

if config.get("print_sample"):
    logger.info(f"Loaded {len(SAMPLES)} samples, {len(TRAITS)} traits")

# ---- Step 9: 包含规则 + TargetPlan 唯一调用点 + rule all ----
include: "rules/01.common.smk"
include: "rules/02.input_gate.smk"
include: "rules/03.genotype_qc.smk"
include: "rules/04.phenotype_qc.smk"
include: "rules/05.population_structure.smk"
include: "rules/06.kinship.smk"
include: "rules/07.association.smk"
include: "rules/08.result_harmonization.smk"
include: "rules/09.locus_definition.smk"
include: "rules/10.local_ld.smk"
include: "rules/11.candidate_genes.smk"
include: "rules/12.visualization.smk"
include: "rules/13.export_summary.smk"
include: "rules/14.deliver.smk"
include: "rules/15.Report.smk"

# TargetPlan 纯函数一次计算（唯一调用点；deliver/report 规则复用计划字段）
ALL_TARGETS = DataDeliver(config=config, samples=SAMPLES, traits=TRAITS)
MODULE_STATUS = ALL_TARGETS.module_status

if config.get("print_target"):
    logger.info("Generated Target Files:\n" + "\n".join(ALL_TARGETS.plan_all))

# Target rule
rule all:
    input:
        ALL_TARGETS.plan_all,
