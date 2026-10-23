#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# GWASFlow - 04.phenotype_qc.smk: 性状级输入准备（每性状独立并行）
#
# 参数全部由 traits.tsv 单点派生（spec §6.3-9，禁止多处重复维护）。
# n_pcs>0 且 population_structure 启用时进 PCA 输入；n_pcs=0 不要求文件存在。

_QC_MIN_SAMPLES = config["parameter"]["qc"].get("min_samples_after_qc", 20)


def _needs_covariates(wildcards):
    trait_def = TRAITS[wildcards.trait]
    return bool(trait_def["numeric_covariates"] or trait_def["categorical_covariates"])


def _needs_pcs_input(wildcards):
    trait_def = TRAITS[wildcards.trait]
    return (
        trait_def["n_pcs"] > 0
        and MODULE_STATUS["population_structure"].enabled
    )


def _needs_skip_pcs(wildcards):
    """only_qc 截断下 n_pcs>0 的性状降级跳过 PC（D3，记录 pc_skipped）。"""
    return (
        TRAITS[wildcards.trait]["n_pcs"] > 0
        and not MODULE_STATUS["population_structure"].enabled
    )


rule prepare_trait:
    """每性状输入准备：FAM 顺序权威对齐 + 全矩阵缺失掩码 + 八件产物。"""
    input:
        gate="01.gate/input_gate.json",
        fam="03.genotype_qc/cohort.fam",
        samples=config["sample_csv"],
        phenotypes=config["phenotypes_tsv"],
        covariates=lambda wc: (
            config["covariates_tsv"] if _needs_covariates(wc) else []
        ),
        pcs=lambda wc: (
            "04.structure/pca.eigenvec" if _needs_pcs_input(wc) else []
        ),
    output:
        pheno="05.phenotype/{trait}/phenotype.gemma.txt",
        covar="05.phenotype/{trait}/covariates.gemma.txt",
        kinship_all="05.phenotype/{trait}/kinship_all.pheno.txt",
        keep="05.phenotype/{trait}/keep.tsv",
        inclusion="05.phenotype/{trait}/sample_inclusion.tsv",
        order="05.phenotype/{trait}/sample_order.tsv",
        design="05.phenotype/{trait}/design_columns.tsv",
        metadata="05.phenotype/{trait}/metadata.json",
    params:
        script=workflow.source_path("../scripts/prepare_trait_inputs.py"),
        trait=lambda wc: wc.trait,
        phenotype_column=lambda wc: TRAITS[wc.trait]["phenotype_column"],
        transform=lambda wc: TRAITS[wc.trait]["transform"],
        unit=lambda wc: TRAITS[wc.trait]["unit"],
        numeric_covariates=lambda wc: ",".join(TRAITS[wc.trait]["numeric_covariates"]),
        categorical_covariates=lambda wc: ",".join(TRAITS[wc.trait]["categorical_covariates"]),
        n_pcs=lambda wc: TRAITS[wc.trait]["n_pcs"],
        skip_pcs=lambda wc: _needs_skip_pcs(wc),
        pcs_path=lambda wc: (
            "04.structure/pca.eigenvec"
            if (TRAITS[wc.trait]["n_pcs"] > 0 and MODULE_STATUS["population_structure"].enabled)
            else ""
        ),
        min_samples=_QC_MIN_SAMPLES,
    log:
        "logs/04.phenotype_qc/{trait}.log",
    benchmark:
        "benchmarks/04.phenotype_qc/{trait}.txt",
    threads: 1,
    resources:
        **rule_resource(config, "low_resource", skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/gwasflow-python.yaml"),
    shell:
        r"""
        set -euo pipefail
        TMP_OUT="05.phenotype/{wildcards.trait}/.tmp.prep.$$"
        rm -rf "$TMP_OUT"; mkdir -p "$(dirname "$TMP_OUT")"
        SKIP_PCS_FLAG=""
        if [ "{params.skip_pcs}" = "True" ]; then SKIP_PCS_FLAG="--skip-pcs"; fi
        PCS_INPUT=""
        if [ -n "{params.pcs_path:q}" ] && [ "{params.skip_pcs}" != "True" ]; then
            PCS_INPUT="--pcs {params.pcs_path:q} --n-pcs {params.n_pcs}"
        fi
        COV_ARG=""
        if [ -n "{input.covariates:q}" ]; then
            COV_ARG="--covariates {input.covariates:q} --numeric-covariates {params.numeric_covariates:q} --categorical-covariates {params.categorical_covariates:q}"
        fi
        python3 {params.script:q} \
            --fam {input.fam:q} \
            --samples {input.samples:q} \
            --phenotypes {input.phenotypes:q} \
            --trait {params.trait:q} \
            --phenotype-column {params.phenotype_column:q} \
            $COV_ARG \
            $PCS_INPUT \
            $SKIP_PCS_FLAG \
            --transform {params.transform:q} \
            --unit {params.unit:q} \
            --min-samples {params.min_samples} \
            --out-dir "$TMP_OUT" > {log:q} 2>&1
        mkdir -p "$(dirname {output.pheno:q})"
        mv "$TMP_OUT/phenotype.gemma.txt" {output.pheno:q}
        mv "$TMP_OUT/covariates.gemma.txt" {output.covar:q}
        mv "$TMP_OUT/kinship_all.pheno.txt" {output.kinship_all:q}
        mv "$TMP_OUT/keep.tsv" {output.keep:q}
        mv "$TMP_OUT/sample_inclusion.tsv" {output.inclusion:q}
        mv "$TMP_OUT/sample_order.tsv" {output.order:q}
        mv "$TMP_OUT/design_columns.tsv" {output.design:q}
        mv "$TMP_OUT/metadata.json" {output.metadata:q}
        rmdir "$TMP_OUT"
        """
