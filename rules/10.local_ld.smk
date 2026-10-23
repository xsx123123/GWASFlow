#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# GWASFlow - 10.local_ld.smk: Regional pairwise LD (per trait; empty loci skip legally)

_LD_PARAMS = config["parameter"]["local_ld"]


rule regional_ld:
    """Regional pairwise LD per locus (on the trait's effective samples; tool filter values audited).

    Input: cohort trio + 05.phenotype/{trait}/keep.tsv + 09.loci/{trait}/loci.tsv.
    Process: for each locus region, run plink2 --r2-unphased (or --r2-phased cols=+dprimeabs for the dprime metric) restricted to the region on the trait's keep set; the tool's own --maf/--geno filter values are audited; per-locus tool failures are recorded in the status file (transparent passthrough, never silently dropped).
    Output: 10.local_ld/{trait}/regional_ld.tsv + regional_ld_status.tsv + ld_index.json.
    """
    input:
        bed="03.genotype_qc/cohort.bed",
        bim="03.genotype_qc/cohort.bim",
        fam="03.genotype_qc/cohort.fam",
        keep="05.phenotype/{trait}/keep.tsv",
        loci="09.loci/{trait}/loci.tsv",
    output:
        ld="10.local_ld/{trait}/regional_ld.tsv",
        status="10.local_ld/{trait}/regional_ld_status.tsv",
        index="10.local_ld/{trait}/ld_index.json",
    params:
        script=workflow.source_path("../scripts/run_regional_ld.py"),
        trait=lambda wc: wc.trait,
        maf=_LD_PARAMS.get("maf", 0.05),
        missing_max=_LD_PARAMS.get("missing_max", 0.10),
        het_max=_LD_PARAMS.get("het_max", 1.0),
        metric=_LD_PARAMS.get("metric", "r2"),
        max_plots=_LD_PARAMS.get("max_regional_plots", 20),
        threads=config["parameter"]["threads"]["plink2"],
    log:
        "logs/10.local_ld/{trait}.log",
    benchmark:
        "benchmarks/10.local_ld/{trait}.txt",
    threads: lambda wc: config["parameter"]["threads"]["plink2"],
    resources:
        **rule_resource(config, "medium_resource", skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/plink2.yaml"),
    shell:
        r"""
        set -euo pipefail
        TMP_OUT="10.local_ld/{wildcards.trait}/.tmp.ld.$$"
        rm -rf "$TMP_OUT"; mkdir -p "$(dirname "$TMP_OUT")"
        python3 {params.script:q} \
            --bfile 03.genotype_qc/cohort \
            --keep {input.keep:q} \
            --loci {input.loci:q} \
            --trait {params.trait:q} \
            --maf {params.maf} \
            --missing-max {params.missing_max} \
            --het-max {params.het_max} \
            --metric {params.metric:q} \
            --max-regional-plots {params.max_plots} \
            --threads {params.threads} \
            --out-dir "$TMP_OUT" > {log:q} 2>&1
        mkdir -p "$(dirname {output.ld:q})"
        mv "$TMP_OUT/regional_ld.tsv" {output.ld:q}
        mv "$TMP_OUT/regional_ld_status.tsv" {output.status:q}
        mv "$TMP_OUT/ld_index.json" {output.index:q}
        rm -rf "$TMP_OUT"
        """
