#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# GWASFlow - 12.visualization.smk: 曼哈顿/QQ/诊断图（每性状 + 项目级）
#
# 契约（spec §5.10）：
# - 图与统一统计/阈值一致（阈值从 metadata 读取）
# - 显著点不因抽样丢失；p_underflow 截断标注不改统计表
# - 上游产物缺失（模块关闭）→ 不生成对应图，plot_manifest 标注原因
# - plot_manifest 与实际图文件一一对应

_PLOT_PARAMS = config["parameter"]["plotting"]
_PLOT_FORMATS = [f.strip() for f in _PLOT_PARAMS.get("formats", "svg,pdf,png").split(",") if f.strip()]
_PLOT_SEED = config.get("seed", 20261023)
_REFERENCE = config.get("GWAS_reference", {}).get(str(config.get("Genome_Version", "")).strip(), {})


_PLOT_GWAS_OUTPUTS = {}
for _fmt in _PLOT_FORMATS:
    _PLOT_GWAS_OUTPUTS[f"manhattan_{_fmt}"] = f"12.plots/{{trait}}/manhattan.{_fmt}"
    _PLOT_GWAS_OUTPUTS[f"qq_{_fmt}"] = f"12.plots/{{trait}}/qq.{_fmt}"
_PLOT_GWAS_OUTPUTS["plot_data"] = "12.plots/{trait}/plot_data.tsv"


rule plot_gwas:
    """曼哈顿/QQ 图（消费统一统计与 metadata 阈值；三格式输出）。"""
    input:
        sumstats="08.results/{trait}/gwas.tsv",
        metadata="08.results/{trait}/metadata.json",
        chromosome_table=(_REFERENCE.get("chromosome_table") or []),
    output:
        **_PLOT_GWAS_OUTPUTS,
    params:
        script=workflow.source_path("../scripts/plot_gwas.R"),
        trait=lambda wc: wc.trait,
        model="gemma_lmm",
        point_budget=_PLOT_PARAMS.get("point_budget", 200000),
        p_zero_display_cap=_PLOT_PARAMS.get("p_zero_display_cap", "1e-300"),
        formats=_PLOT_PARAMS.get("formats", "svg,pdf,png"),
        seed=_PLOT_SEED,
    log:
        "logs/12.visualization/{trait}_gwas.log",
    benchmark:
        "benchmarks/12.visualization/{trait}_gwas.txt",
    threads: lambda wc: config["parameter"]["threads"]["plotting"],
    resources:
        **rule_resource(config, "low_resource", skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/r-plot.yaml"),
    shell:
        r"""
        set -euo pipefail
        mkdir -p 12.plots/{wildcards.trait}
        CT_ARG=""
        if [ -n "{input.chromosome_table:q}" ]; then CT_ARG="--chromosome-table {input.chromosome_table:q}"; fi
        Rscript {params.script:q} \
            --sumstats {input.sumstats:q} \
            --metadata {input.metadata:q} \
            $CT_ARG \
            --trait {params.trait:q} \
            --model {params.model:q} \
            --point-budget {params.point_budget} \
            --preserve-significant true \
            --p-zero-display-cap {params.p_zero_display_cap:q} \
            --formats {params.formats:q} \
            --seed {params.seed} \
            --out-dir 12.plots/{wildcards.trait} > {log:q} 2>&1
        """


rule plot_trait_diag:
    """每性状诊断图：表型分布 + 样本纳入概况（消费真实准备产物）。"""
    input:
        phenotype="05.phenotype/{trait}/phenotype.gemma.txt",
        inclusion="05.phenotype/{trait}/sample_inclusion.tsv",
    output:
        pheno_dist="12.plots/{trait}/phenotype_distribution.png",
        inclusion_plot="12.plots/{trait}/sample_inclusion.png",
    params:
        script=workflow.source_path("../scripts/plot_diagnostics.R"),
        trait=lambda wc: wc.trait,
        unit=lambda wc: TRAITS[wc.trait]["unit"],
    log:
        "logs/12.visualization/{trait}_diag.log",
    benchmark:
        "benchmarks/12.visualization/{trait}_diag.txt",
    threads: lambda wc: config["parameter"]["threads"]["plotting"],
    resources:
        **rule_resource(config, "low_resource", skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/r-plot.yaml"),
    shell:
        r"""
        set -euo pipefail
        mkdir -p 12.plots/{wildcards.trait}
        Rscript {params.script:q} \
            --trait {params.trait:q} \
            --unit {params.unit:q} \
            --phenotype {input.phenotype:q} \
            --sample-inclusion {input.inclusion:q} \
            --out-dir 12.plots/{wildcards.trait} > {log:q} 2>&1
        """


def _project_diag_needed():
    return (
        MODULE_STATUS["population_structure"].enabled
        and MODULE_STATUS["kinship_analysis"].enabled
    )


_K_MATRIX_NAME = "kinship.cXX.txt" if config["parameter"]["kinship"].get("type", "centered") == "centered" else "kinship.sXX.txt"


rule plot_project_diag:
    """项目级诊断图：PCA + K 热图（仅结构与 kinship 模块均启用时）。"""
    input:
        pcs="04.structure/pca.eigenvec",
        kinship=f"06.kinship/{_K_MATRIX_NAME}",
    output:
        pca="12.plots/diagnostics/pca.png",
        kinship_heatmap="12.plots/diagnostics/kinship_heatmap.png",
    params:
        script=workflow.source_path("../scripts/plot_diagnostics.R"),
    log:
        "logs/12.visualization/project_diag.log",
    benchmark:
        "benchmarks/12.visualization/project_diag.txt",
    threads: lambda wc: config["parameter"]["threads"]["plotting"],
    resources:
        **rule_resource(config, "low_resource", skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/r-plot.yaml"),
    shell:
        r"""
        set -euo pipefail
        mkdir -p 12.plots/diagnostics
        Rscript {params.script:q} \
            --trait project \
            --pcs {input.pcs:q} \
            --kinship {input.kinship:q} \
            --out-dir 12.plots/diagnostics > {log:q} 2>&1
        """


def _plot_manifest_inputs(wildcards):
    inputs = []
    for fmt in _PLOT_FORMATS:
        inputs.extend(f"12.plots/{t}/manhattan.{fmt}" for t in TRAITS.keys())
        inputs.extend(f"12.plots/{t}/qq.{fmt}" for t in TRAITS.keys())
    inputs.extend(f"12.plots/{t}/phenotype_distribution.png" for t in TRAITS.keys())
    inputs.extend(f"12.plots/{t}/sample_inclusion.png" for t in TRAITS.keys())
    inputs.extend(f"12.plots/{t}/plot_data.tsv" for t in TRAITS.keys())
    if _project_diag_needed():
        inputs.extend(["12.plots/diagnostics/pca.png",
                       "12.plots/diagnostics/kinship_heatmap.png"])
    return inputs


def _plot_manifest_not_generated():
    reasons = []
    if not _project_diag_needed():
        reasons.append(
            "12.plots/diagnostics/pca.png：population_structure 或 kinship_analysis "
            "模块未启用（上游产物不存在，不生成并如实标注）"
        )
        reasons.append(
            "12.plots/diagnostics/kinship_heatmap.png：同上"
        )
    return reasons


rule plot_manifest:
    """图表清单：图身份、来源数据、格式（与实际图文件一一对应）。"""
    input:
        _plot_manifest_inputs,
    output:
        manifest="12.plots/plot_manifest.json",
    params:
        script=workflow.source_path("../scripts/build_plot_manifest.py"),
        plot_files=lambda wc: _plot_manifest_inputs(None),
        not_generated_b64=lambda wc: __import__("base64").b64encode(
            "\n".join(_plot_manifest_not_generated()).encode()
        ).decode(),
    log:
        "logs/12.visualization/plot_manifest.log",
    benchmark:
        "benchmarks/12.visualization/plot_manifest.txt",
    threads: 1,
    resources:
        **rule_resource(config, "low_resource", skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/gwasflow-python.yaml"),
    shell:
        r"""
        set -euo pipefail
        TMP_DIR=".tmp.plotmanifest.$$"
        rm -rf "$TMP_DIR"; mkdir -p "$TMP_DIR"
        echo {params.not_generated_b64:q} | base64 -d > "$TMP_DIR/not_generated.txt"
        python3 {params.script:q}             --not-generated-file "$TMP_DIR/not_generated.txt"             $(for f in {params.plot_files}; do printf ' --plot-file %s' "$f"; done)             --out {output.manifest:q} > {log:q} 2>&1
        rm -rf "$TMP_DIR"
        """
