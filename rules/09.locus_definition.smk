#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# GWASFlow - 09.locus_definition.smk: Significant locus definition (per trait)

_LOCI_PARAMS = config["parameter"]["loci"]
_ASSOC_PARAMS = config["parameter"]["association"]


def _assoc_test_file(wildcards):
    test = _ASSOC_PARAMS.get("primary_test", "wald")
    mapping = {
        "wald": "association.assoc.txt",
        "lrt": "association.lrt.txt",
        "score": "association.score.txt",
        "all": "association.assoc.txt",
    }
    return f"08.results/{{trait}}/metadata.json"


rule define_loci:
    """Locus definition: lead selection (threshold from metadata as the only source) + clump de-redundancy + region boundary registration.

    Input: 08.results/{trait}/gwas.tsv + significant.tsv + metadata.json (threshold) + cohort trio + trait keep.tsv.
    Process: select lead variants from the significant set using the metadata threshold; collapse redundant signals via plink2 clump (r2/kb); register locus boundaries (LD-based with fixed fallback flank); an empty significant set legally yields empty loci.
    Output: 09.loci/{trait}/lead.tsv + loci.tsv + locus_manifest.json.
    """
    input:
        sumstats="08.results/{trait}/gwas.tsv",
        significant="08.results/{trait}/significant.tsv",
        metadata="08.results/{trait}/metadata.json",
        bed="03.genotype_qc/cohort.bed",
        bim="03.genotype_qc/cohort.bim",
        fam="03.genotype_qc/cohort.fam",
        keep="05.phenotype/{trait}/keep.tsv",
    output:
        lead="09.loci/{trait}/lead.tsv",
        loci="09.loci/{trait}/loci.tsv",
        manifest="09.loci/{trait}/locus_manifest.json",
    params:
        script=workflow.source_path("../scripts/define_loci.py"),
        trait=lambda wc: wc.trait,
        model="gemma_lmm",
        clump_r2=_LOCI_PARAMS.get("clump_r2", 0.20),
        clump_kb=_LOCI_PARAMS.get("clump_kb", 500),
        region_boundary=_LOCI_PARAMS.get("region_boundary", "ld"),
        fallback_flank_bp=_LOCI_PARAMS.get("fallback_flank_bp", 50000),
        threads=config["parameter"]["threads"]["plink2"],
    log:
        "logs/09.locus_definition/{trait}.log",
    benchmark:
        "benchmarks/09.locus_definition/{trait}.txt",
    threads: lambda wc: config["parameter"]["threads"]["plink2"],
    resources:
        **rule_resource(config, "medium_resource", skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/plink2.yaml"),
    shell:
        r"""
        set -euo pipefail
        TMP_OUT="09.loci/{wildcards.trait}/.tmp.loci.$$"
        rm -rf "$TMP_OUT"; mkdir -p "$(dirname "$TMP_OUT")"
        python3 {params.script:q} \
            --sumstats {input.sumstats:q} \
            --bfile 03.genotype_qc/cohort \
            --significant {input.significant:q} \
            --association-metadata {input.metadata:q} \
            --trait {params.trait:q} \
            --model {params.model:q} \
            --keep {input.keep:q} \
            --mode ld_clump \
            --clump-r2 {params.clump_r2} \
            --clump-kb {params.clump_kb} \
            --region-boundary {params.region_boundary:q} \
            --fallback-flank-bp {params.fallback_flank_bp} \
            --threads {params.threads} \
            --out-dir "$TMP_OUT" > {log:q} 2>&1
        mkdir -p "$(dirname {output.lead:q})"
        mv "$TMP_OUT/lead.tsv" {output.lead:q}
        mv "$TMP_OUT/loci.tsv" {output.loci:q}
        mv "$TMP_OUT/locus_manifest.json" {output.manifest:q}
        rm -rf "$TMP_OUT"
        """
