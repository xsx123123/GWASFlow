#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# GWASFlow - 03.genotype_qc.smk: 分阶段基因型 QC（四阶段独立执行）
#
# cohort 三件套的 FAM 顺序成为后续全流程唯一对齐权威。

_QC_PARAMS = config["parameter"]["qc"]


rule genotype_qc:
    """四阶段 QC：pre_geno → mind → geno/maf/mac → final_stats（独立 PLINK2 调用）。"""
    input:
        bed="02.genotype/raw_engine.bed",
        bim="02.genotype/raw_engine.bim",
        fam="02.genotype/raw_engine.fam",
        gate="01.gate/input_gate.json",
    output:
        bed="03.genotype_qc/cohort.bed",
        bim="03.genotype_qc/cohort.bim",
        fam="03.genotype_qc/cohort.fam",
        stage1_stats="03.genotype_qc/stage1.stats.tsv",
        stage2_stats="03.genotype_qc/stage2.stats.tsv",
        stage3_stats="03.genotype_qc/stage3.stats.tsv",
        removed_v1="03.genotype_qc/removed_variants_stage1.tsv",
        removed_s2="03.genotype_qc/removed_samples_stage2.tsv",
        removed_v3="03.genotype_qc/removed_variants_stage3.tsv",
        final_stats_tsv="03.genotype_qc/final_stats.tsv",
        final_stats_json="03.genotype_qc/final_stats.json",
        audit="03.genotype_qc/qc_audit.json",
    params:
        script=workflow.source_path("../scripts/genotype_qc.py"),
        chr_count=_QC_PARAMS.get("chr_count", 9),
        pre_geno=_QC_PARAMS.get("pre_geno", 0.20),
        mind=_QC_PARAMS.get("mind", 0.10),
        geno=_QC_PARAMS.get("geno", 0.10),
        maf=_QC_PARAMS.get("maf", 0.05),
        mac=_QC_PARAMS.get("mac", 5),
        min_samples=_QC_PARAMS.get("min_samples_after_qc", 20),
        threads=config["parameter"]["threads"]["plink2"],
        memory_mb=_QC_PARAMS.get("memory_mb", 8000),
    log:
        "logs/03.genotype_qc/genotype_qc.log",
    benchmark:
        "benchmarks/03.genotype_qc/genotype_qc.txt",
    threads: lambda wc: config["parameter"]["threads"]["plink2"],
    resources:
        **rule_resource(config, "medium_resource", skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/plink2.yaml"),
    shell:
        r"""
        set -euo pipefail
        TMP_OUT="03.genotype_qc/.tmp.qc.$$"
        rm -rf "$TMP_OUT"; mkdir -p "$(dirname "$TMP_OUT")"
        python3 {params.script:q} \
            --mode qc \
            --bfile 02.genotype/raw_engine \
            --gate {input.gate:q} \
            --chr-count {params.chr_count} \
            --pre-geno {params.pre_geno} \
            --mind {params.mind} \
            --geno {params.geno} \
            --maf {params.maf} \
            --mac {params.mac} \
            --min-samples-after-qc {params.min_samples} \
            --threads {params.threads} \
            --memory-mb {params.memory_mb} \
            --out-dir "$TMP_OUT" > {log:q} 2>&1
        mkdir -p "$(dirname {output.bed:q})"
        mv "$TMP_OUT/cohort.bed" {output.bed:q}
        mv "$TMP_OUT/cohort.bim" {output.bim:q}
        mv "$TMP_OUT/cohort.fam" {output.fam:q}
        mv "$TMP_OUT/stage1.stats.tsv" {output.stage1_stats:q}
        mv "$TMP_OUT/stage2.stats.tsv" {output.stage2_stats:q}
        mv "$TMP_OUT/stage3.stats.tsv" {output.stage3_stats:q}
        mv "$TMP_OUT/removed_variants_stage1.tsv" {output.removed_v1:q}
        mv "$TMP_OUT/removed_samples_stage2.tsv" {output.removed_s2:q}
        mv "$TMP_OUT/removed_variants_stage3.tsv" {output.removed_v3:q}
        mv "$TMP_OUT/final_stats.tsv" {output.final_stats_tsv:q}
        mv "$TMP_OUT/final_stats.json" {output.final_stats_json:q}
        mv "$TMP_OUT/qc_audit.json" {output.audit:q}
        rm -rf "$TMP_OUT"
        """
