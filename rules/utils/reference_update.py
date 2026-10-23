#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GWASFlow 参考路径解析与按需校验。

resolve_reference_paths: 将 GWAS_reference 中相对 reference_path 的路径
拼接为绝对路径（原地更新，仅字符串处理不做 IO）。
check_gwas_reference_paths: 按生效模块子集校验所需参考资源存在性
（仅启用模块所需；only_qc 且未启用候选基因时 GFF3 缺失不阻断）。
"""

import os
import sys


def _die(message: str, code: int = 2) -> None:
    print(f"[reference_update ERROR] {message}", file=sys.stderr)
    sys.exit(code)


def resolve_reference_paths(config, index_keys=None, base_path=None, logger=None):
    """解析 index_keys（默认 ["GWAS_reference"]）下各版本的相对路径为绝对路径。

    契约对齐 FlowFrame 模板：value 不以根目录开头且非绝对路径时拼接。
    """
    if index_keys is None:
        index_keys = ["GWAS_reference"]

    root_dir = base_path or config.get("reference_path")
    if not root_dir:
        return

    for section_key in index_keys:
        section = config.get(section_key, {})
        if not isinstance(section, dict):
            continue
        for version, ref_dict in section.items():
            if not isinstance(ref_dict, dict):
                continue
            for key, value in ref_dict.items():
                if isinstance(value, str):
                    if not value.startswith(root_dir) and not value.startswith("/"):
                        ref_dict[key] = os.path.join(root_dir, value)


# 各模块所需参考资源键（按需校验子集依据）
_MODULE_REFERENCE_REQUIREMENTS = {
    "input_gate": [],  # 门控只登记，不强校验具体资源
    "genotype_normalize": ["genome_fasta", "genome_fai", "chromosome_rename_map"],
    "genotype_qc": [],
    "population_structure": [],
    "phenotype_qc": [],
    "kinship_analysis": [],
    "association_analysis": [],
    "result_harmonization": [],
    "locus_analysis": [],
    "local_ld": [],
    "candidate_genes": ["gff3"],
    "visualization": ["chromosome_table"],
    "export_summary": [],
}


def check_gwas_reference_paths(config, enabled_modules, logger=None):
    """按生效模块子集校验当前 Genome_Version 的参考资源。

    only_qc 截断后未启用的模块其资源缺失不阻断（spec §5.1.1-4）。
    缺失时列出 genome/key/path 明细并退出 1。
    """
    version = str(config.get("Genome_Version", "")).strip()
    registry = config.get("GWAS_reference", {})
    ref_dict = registry.get(version)
    if not isinstance(ref_dict, dict) or not ref_dict:
        _die(
            f"GWAS_reference 中未登记当前 Genome_Version '{version}' 的参考资源"
        )

    required_keys = set()
    for module in enabled_modules:
        for key in _MODULE_REFERENCE_REQUIREMENTS.get(module, []):
            required_keys.add(key)

    if "genotype_normalize" in enabled_modules:
        fmt = config.get("genotype_format")
        if fmt in ("vcf", "bcf"):
            required_keys.update(["genome_fasta", "genome_fai", "chromosome_rename_map"])
        else:
            # BED 路线：BIM 已是引擎数值染色体，无需 FASTA/FAI 与重命名映射
            required_keys.discard("genome_fasta")
            required_keys.discard("genome_fai")
            required_keys.discard("chromosome_rename_map")

    missing = []
    for key in sorted(required_keys):
        path = ref_dict.get(key)
        if not path or not os.path.exists(path):
            missing.append((version, key, path or "<未配置>"))

    if missing:
        lines = [f"  - {g} / {k}: {p} 不存在" for g, k, p in missing]
        msg = "参考资源缺失（按生效模块所需）:\n" + "\n".join(lines)
        if logger:
            logger.error(msg)
        _die(msg, code=1)

    if logger:
        logger.info(
            f"参考资源校验通过: {version}（所需键 {sorted(required_keys)}）"
        )
