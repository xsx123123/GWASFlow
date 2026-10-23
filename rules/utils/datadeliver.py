#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GWASFlow 模块目标收集器（trait 维度展开）。

硬约束（框架规范 §5.6）：
- 每条返回路径与对应 rule output 逐字一致（改 output 必须同 PR 改收集器）
- 空壳收集器视同缺陷
- 运行后数量型目标（locus 数量未知）按每性状批量区域目录进 manifest，
  禁止在纯目标函数中读取尚不存在的统计表（D7）
- 收集器只读 config/samples/traits/module_status（纯函数）
"""

from snakemake.io import expand


def gate_outputs(config, samples, traits, module_status):
    """01.gate/ 三件 + effective_config.yaml（02 verify_inputs 规则 output）。"""
    return [
        "01.gate/input_gate.json",
        "01.gate/md5_check.tsv",
        "01.gate/reference_check.json",
        "13.summary/effective_config.yaml",
    ]


def normalize_outputs(config, samples, traits, module_status):
    """02.genotype/ 规范化产物（raw_engine 三件套 + variant_map + 审计）。"""
    return [
        "02.genotype/raw_engine.bed",
        "02.genotype/raw_engine.bim",
        "02.genotype/raw_engine.fam",
        "02.genotype/variant_map.tsv",
        "02.genotype/normalize_audit.json",
    ]


def qc_outputs(config, samples, traits, module_status):
    """03.genotype_qc/ 四阶段 QC 全套产物。"""
    return [
        "03.genotype_qc/cohort.bed",
        "03.genotype_qc/cohort.bim",
        "03.genotype_qc/cohort.fam",
        "03.genotype_qc/stage1.stats.tsv",
        "03.genotype_qc/stage2.stats.tsv",
        "03.genotype_qc/stage3.stats.tsv",
        "03.genotype_qc/removed_variants_stage1.tsv",
        "03.genotype_qc/removed_samples_stage2.tsv",
        "03.genotype_qc/removed_variants_stage3.tsv",
        "03.genotype_qc/final_stats.tsv",
        "03.genotype_qc/final_stats.json",
        "03.genotype_qc/qc_audit.json",
    ]


def structure_outputs(config, samples, traits, module_status):
    """04.structure/ 剪枝集与 PCA 产物。"""
    return [
        "04.structure/pca_prune.prune.in",
        "04.structure/pca_prune.prune.out",
        "04.structure/pca.eigenvec",
        "04.structure/pca.eigenval",
    ]


TRAIT_PREP_FILES = [
    "phenotype.gemma.txt",
    "covariates.gemma.txt",
    "kinship_all.pheno.txt",
    "keep.tsv",
    "sample_inclusion.tsv",
    "sample_order.tsv",
    "design_columns.tsv",
    "metadata.json",
]


def phenotype_outputs(config, samples, traits, module_status):
    """05.phenotype/{trait}/ 每性状八件产物。"""
    if not traits:
        return []
    return expand("05.phenotype/{trait}/" + "{file}", trait=list(traits.keys()), file=TRAIT_PREP_FILES)


def kinship_outputs(config, samples, traits, module_status):
    """06.kinship/ 全队列 K + 顺序 + 验证报告（文件名按 K 类型）。"""
    k_type = (config.get("parameter", {}).get("kinship", {}).get("type") or "centered")
    matrix_name = "kinship.cXX.txt" if k_type == "centered" else "kinship.sXX.txt"
    return [
        f"06.kinship/{matrix_name}",
        "06.kinship/kinship.log.txt",
        "06.kinship/sample_order.tsv",
        "06.kinship/kinship_validation.json",
    ]


def association_outputs(config, samples, traits, module_status):
    """07.association/{trait}/ 每性状 MAC + GEMMA 原始输出 + run_metadata。

    检验输出文件名按 primary_test 展开（wald=assoc/lrt/score；all=三件独立归档）。
    """
    if not traits:
        return []
    test = (config.get("parameter", {}).get("association", {})
            .get("primary_test", "wald"))
    assoc_files = {
        "wald": ["association.assoc.txt"],
        "lrt": ["association.lrt.txt"],
        "score": ["association.score.txt"],
        "all": ["association.assoc.txt", "association.lrt.txt",
                "association.score.txt"],
    }[test]
    files = ["trait_mac.tsv", "trait_mac_summary.json"] + assoc_files + [
        "association.log.txt", "gemma.help.txt", "run_metadata.json",
    ]
    return expand("07.association/{trait}/" + "{file}", trait=list(traits.keys()), file=files)


def harmonization_outputs(config, samples, traits, module_status):
    """08.results/{trait}/ 每性状统一统计五件套。"""
    if not traits:
        return []
    files = [
        "gwas.tsv",
        "significant.tsv",
        "invalid_variants.tsv",
        "unreported_variants.tsv",
        "metadata.json",
    ]
    return expand("08.results/{trait}/" + "{file}", trait=list(traits.keys()), file=files)


def locus_outputs(config, samples, traits, module_status):
    """09.loci/{trait}/ 每性状 lead/loci/manifest。"""
    if not traits:
        return []
    files = ["lead.tsv", "loci.tsv", "locus_manifest.json"]
    return expand("09.loci/{trait}/" + "{file}", trait=list(traits.keys()), file=files)


def local_ld_outputs(config, samples, traits, module_status):
    """10.local_ld/{trait}/ 每性状区域 LD 数据/状态/索引。"""
    if not traits:
        return []
    files = ["regional_ld.tsv", "regional_ld_status.tsv", "ld_index.json"]
    return expand("10.local_ld/{trait}/" + "{file}", trait=list(traits.keys()), file=files)


def candidate_outputs(config, samples, traits, module_status):
    """11.candidates/{trait}/ 每性状候选基因表与元数据。"""
    if not traits:
        return []
    files = ["candidates.tsv", "candidates.tsv.meta.json"]
    return expand("11.candidates/{trait}/" + "{file}", trait=list(traits.keys()), file=files)


PLOT_FORMATS = ["svg", "pdf", "png"]


def visualization_outputs(config, samples, traits, module_status):
    """12.plots/ 每性状曼哈顿/QQ/表型诊断 + 项目级诊断与 plot_manifest。"""
    outputs = []
    if traits:
        for fmt in PLOT_FORMATS:
            outputs.extend(
                expand(f"12.plots/{{trait}}/manhattan.{fmt}", trait=list(traits.keys()))
            )
            outputs.extend(
                expand(f"12.plots/{{trait}}/qq.{fmt}", trait=list(traits.keys()))
            )
        outputs.extend(
            expand("12.plots/{trait}/phenotype_distribution.png", trait=list(traits.keys()))
        )
        outputs.extend(
            expand("12.plots/{trait}/sample_inclusion.png", trait=list(traits.keys()))
        )
        outputs.extend(
            expand("12.plots/{trait}/plot_data.tsv", trait=list(traits.keys()))
        )
    if (
        module_status.get("population_structure") and module_status["population_structure"].enabled
        and module_status.get("kinship_analysis") and module_status["kinship_analysis"].enabled
    ):
        outputs.extend(["12.plots/diagnostics/pca.png", "12.plots/diagnostics/kinship_heatmap.png"])
    outputs.append("12.plots/plot_manifest.json")
    return outputs


def summary_outputs(config, samples, traits, module_status):
    """13.summary/ 项目摘要与显式产物索引。"""
    return [
        "13.summary/project_summary.json",
        "13.summary/artifact_index.tsv",
    ]


def delivery_outputs(config, samples, traits, module_status):
    """14.data_deliver/ 交付 manifest 三件套。"""
    return [
        "14.data_deliver/delivery_manifest.json",
        "14.data_deliver/delivery_manifest.md5",
        "14.data_deliver/delivery_details.log",
    ]


def report_outputs(config, samples, traits, module_status):
    """15.report/ 报告产物（report=true 时；渲染器未验收时运行前报错）。"""
    if not module_status.get("report") or not module_status["report"].enabled:
        return []
    return ["15.report/report_manifest.json"]


# 分析类模块（12 个）与收尾模块（summary/deliver/report）
ANALYSIS_COLLECTORS = [
    ("input_gate", gate_outputs),
    ("genotype_normalize", normalize_outputs),
    ("genotype_qc", qc_outputs),
    ("population_structure", structure_outputs),
    ("phenotype_qc", phenotype_outputs),
    ("kinship_analysis", kinship_outputs),
    ("association_analysis", association_outputs),
    ("result_harmonization", harmonization_outputs),
    ("locus_analysis", locus_outputs),
    ("local_ld", local_ld_outputs),
    ("candidate_genes", candidate_outputs),
    ("visualization", visualization_outputs),
]


def collect_all_targets(config, samples, traits, module_status):
    """计算五类目标（analysis/summary/report/delivery → 由调用方组装 all）。

    返回 (analysis_list, summary_list, report_list, delivery_list)。
    """
    analysis = []
    for module, collector in ANALYSIS_COLLECTORS:
        if module_status.get(module) and module_status[module].enabled:
            analysis.extend(collector(config, samples, traits, module_status))

    # summary 类：effective_config（gate 规则 output）+ project_summary + artifact_index
    summary = []
    if module_status.get("input_gate") and module_status["input_gate"].enabled:
        summary.append("13.summary/effective_config.yaml")
    if module_status.get("export_summary") and module_status["export_summary"].enabled:
        summary.extend(summary_outputs(config, samples, traits, module_status))

    report = report_outputs(config, samples, traits, module_status)
    delivery = []
    if module_status.get("deliver") and module_status["deliver"].enabled:
        delivery.extend(delivery_outputs(config, samples, traits, module_status))

    return analysis, summary, report, delivery
