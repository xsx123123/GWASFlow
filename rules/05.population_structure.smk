#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# GWASFlow - 05.population_structure.smk: LD 剪枝 + 剪枝集 PCA
#
# 用途隔离（spec §5.4.1-1）：PCA 使用剪枝标记集；关联分析使用完整 cohort 标记集。
# 规则编号 ≠ 执行顺序：04.phenotype_qc 依赖本规则产物（n_pcs>0 时）。

_QC_PARAMS = config["parameter"]["qc"]
_STRUCT_PARAMS = config["parameter"]["structure"]


rule prune_pca:
    """LD 剪枝（--indep-pairwise 窗口 kb/步长/r2）+ 剪枝集 PCA（scols=fid 保留 FID）。"""
    input:
        bed="03.genotype_qc/cohort.bed",
        bim="03.genotype_qc/cohort.bim",
        fam="03.genotype_qc/cohort.fam",
    output:
        prune_in="04.structure/pca_prune.prune.in",
        prune_out="04.structure/pca_prune.prune.out",
        eigenvec="04.structure/pca.eigenvec",
        eigenval="04.structure/pca.eigenval",
    params:
        script=workflow.source_path("../scripts/genotype_qc.py"),
        chr_count=_QC_PARAMS.get("chr_count", 9),
        prune_window_kb=_STRUCT_PARAMS.get("prune_window_kb", 50),
        prune_step=_STRUCT_PARAMS.get("prune_step", 5),
        prune_r2=_STRUCT_PARAMS.get("prune_r2", 0.20),
        pca_components=_STRUCT_PARAMS.get("pca_components", 10),
        threads=config["parameter"]["threads"]["plink2"],
        memory_mb=_QC_PARAMS.get("memory_mb", 8000),
    log:
        "logs/05.population_structure/prune_pca.log",
    benchmark:
        "benchmarks/05.population_structure/prune_pca.txt",
    threads: lambda wc: config["parameter"]["threads"]["plink2"],
    resources:
        **rule_resource(config, "medium_resource", skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/plink2.yaml"),
    shell:
        r"""
        set -euo pipefail
        TMP_OUT="04.structure/.tmp.structure.$$"
        rm -rf "$TMP_OUT"; mkdir -p "$(dirname "$TMP_OUT")"
        python3 {params.script:q} \
            --mode structure \
            --bfile 03.genotype_qc/cohort \
            --chr-count {params.chr_count} \
            --prune-window-kb {params.prune_window_kb} \
            --prune-step {params.prune_step} \
            --prune-r2 {params.prune_r2} \
            --pca-components {params.pca_components} \
            --threads {params.threads} \
            --memory-mb {params.memory_mb} \
            --out-dir "$TMP_OUT" > {log:q} 2>&1
        mkdir -p "$(dirname {output.prune_in:q})"
        mv "$TMP_OUT/pca_prune.prune.in" {output.prune_in:q}
        mv "$TMP_OUT/pca_prune.prune.out" {output.prune_out:q}
        mv "$TMP_OUT/pca.eigenvec" {output.eigenvec:q}
        mv "$TMP_OUT/pca.eigenval" {output.eigenval:q}
        rm -rf "$TMP_OUT"
        """
