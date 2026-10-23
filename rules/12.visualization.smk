#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# GWASFlow - 12.visualization.smk: Manhattan/QQ/diagnostic plots (per trait + project level)
#
# Contracts (spec §5.10):
# - Plots agree with the harmonized statistics/threshold (threshold read from metadata)
# - Significant points are never lost to sampling; the p_underflow display cap never modifies the statistics table
# - Missing upstream artifacts (module disabled) -> the corresponding plot is not generated and plot_manifest records why
# - plot_manifest maps one-to-one to the actual plot files

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
    """Static Manhattan/QQ plots (consume harmonized statistics and the metadata threshold; three formats).

    Input: 08.results/{trait}/gwas.tsv + metadata.json (+ optional chromosome table).
    Process: render the Manhattan plot (genome-offset x axis from chromosome-table lengths, Bonferroni line from metadata, dotted p_underflow display-cap line when present) and the QQ plot (expected -log10(rank/(n+1)) vs observed); non-significant points beyond point_budget are sampled while significant points are always kept; the sampling is audited in plot_data.tsv.
    Output: 12.plots/{trait}/manhattan.{svg,pdf,png} + qq.{svg,pdf,png} + plot_data.tsv.
    """
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


rule plot_gwas_interactive:
    """Interactive Manhattan/QQ plots (plotly WebGL HTML; for exploration — the static three-format plots remain produced by plot_gwas).

    Input: 08.results/{trait}/gwas.tsv + metadata.json (+ optional chromosome table).
    Process: render per-chromosome WebGL scatter traces via plotly go.Scattergl following the interaction model of the plotly.com Manhattan-plot tutorial (zoom/pan/box/lasso, hover with variant identity, P, beta and effect allele); Bonferroni threshold line from metadata; p_underflow display cap consistent with the static plot; full rendering by default with sampling only above max_render_points (significant points always kept, seed reproducible).
    Output: 12.plots/{trait}/manhattan_interactive.html + qq_interactive.html (offline-openable, plotly.js embedded).
    """
    input:
        sumstats="08.results/{trait}/gwas.tsv",
        metadata="08.results/{trait}/metadata.json",
        chromosome_table=(_REFERENCE.get("chromosome_table") or []),
    output:
        manhattan_html="12.plots/{trait}/manhattan_interactive.html",
        qq_html="12.plots/{trait}/qq_interactive.html",
    params:
        script=workflow.source_path("../scripts/plot_gwas_interactive.py"),
        trait=lambda wc: wc.trait,
        model="gemma_lmm",
        max_render_points=_PLOT_PARAMS.get("max_render_points", 2000000),
        p_zero_display_cap=_PLOT_PARAMS.get("p_zero_display_cap", "1e-300"),
        plotlyjs=_PLOT_PARAMS.get("plotlyjs", "embed"),
        seed=_PLOT_SEED,
    log:
        "logs/12.visualization/{trait}_gwas_interactive.log",
    benchmark:
        "benchmarks/12.visualization/{trait}_gwas_interactive.txt",
    threads: lambda wc: config["parameter"]["threads"]["plotting"],
    resources:
        **rule_resource(config, "low_resource", skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/gwasflow-plotly.yaml"),
    shell:
        r"""
        set -euo pipefail
        mkdir -p 12.plots/{wildcards.trait}
        CT_ARG=""
        if [ -n "{input.chromosome_table:q}" ]; then CT_ARG="--chromosome-table {input.chromosome_table:q}"; fi
        python3 {params.script:q} \
            --sumstats {input.sumstats:q} \
            --metadata {input.metadata:q} \
            $CT_ARG \
            --trait {params.trait:q} \
            --model {params.model:q} \
            --max-render-points {params.max_render_points} \
            --p-zero-display-cap {params.p_zero_display_cap:q} \
            --plotlyjs {params.plotlyjs:q} \
            --seed {params.seed} \
            --out-dir 12.plots/{wildcards.trait} > {log:q} 2>&1
        """


rule plot_trait_diag:
    """Per-trait diagnostic plots: phenotype distribution + sample inclusion overview (consume real preparation artifacts).

    Input: 05.phenotype/{trait}/phenotype.gemma.txt + sample_inclusion.tsv.
    Process: render the phenotype distribution and the sample inclusion summary from the prepared trait inputs.
    Output: 12.plots/{trait}/phenotype_distribution.png + sample_inclusion.png.
    """
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
    """Project-level diagnostic plots: PCA + kinship heatmap (only when both structure and kinship modules are enabled).

    Input: 04.structure/pca.eigenvec + 06.kinship kinship matrix.
    Process: render the PC1/PC2 scatter from the LD-pruned PCA and the kinship heatmap (first 200 samples).
    Output: 12.plots/diagnostics/pca.png + kinship_heatmap.png.
    """
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
    inputs.extend(f"12.plots/{t}/manhattan_interactive.html" for t in TRAITS.keys())
    inputs.extend(f"12.plots/{t}/qq_interactive.html" for t in TRAITS.keys())
    if _project_diag_needed():
        inputs.extend(["12.plots/diagnostics/pca.png",
                       "12.plots/diagnostics/kinship_heatmap.png"])
    return inputs


def _plot_manifest_not_generated():
    reasons = []
    if not _project_diag_needed():
        reasons.append(
            "12.plots/diagnostics/pca.png: population_structure or kinship_analysis "
            "module not enabled (upstream artifacts absent; not generated, recorded honestly)"
        )
        reasons.append(
            "12.plots/diagnostics/kinship_heatmap.png: same as above"
        )
    return reasons


rule plot_manifest:
    """Plot manifest: plot identity, source data, format (one-to-one with the actual plot files).

    Input: all plot files listed by _plot_manifest_inputs (Manhattan/QQ per format, diagnostics, plot data, interactive HTML) + not-generated reasons.
    Process: register each existing plot file with role, format and size; declared-but-missing files abort the manifest (no silent omissions); module-disabled plots are recorded with their reasons.
    Output: 12.plots/plot_manifest.json.
    """
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
