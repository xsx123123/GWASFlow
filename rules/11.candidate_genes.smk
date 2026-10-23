#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# GWASFlow - 11.candidate_genes.smk: Candidate gene annotation (per trait)

_CAND_PARAMS = config["parameter"]["candidate"]
_REFERENCE = config.get("GWAS_reference", {}).get(str(config.get("Genome_Version", "")).strip(), {})


rule annotate_genes:
    """Candidate gene annotation: GFF3 interval intersection + function table join (multi-row genes aggregated first).

    Input: 09.loci/{trait}/loci.tsv + reference GFF3 (+ optional function table).
    Process: intersect locus regions with gene features from the GFF3 (with configurable flank); join gene descriptions from the function table, aggregating multi-row entries before the join.
    Output: 11.candidates/{trait}/candidates.tsv + candidates.tsv.meta.json.
    """
    input:
        loci="09.loci/{trait}/loci.tsv",
        gff3=(_REFERENCE.get("gff3") or []),
        function_table=(_REFERENCE.get("function_table") or []),
    output:
        candidates="11.candidates/{trait}/candidates.tsv",
        meta="11.candidates/{trait}/candidates.tsv.meta.json",
    params:
        script=workflow.source_path("../scripts/annotate_loci_genes.py"),
        genome_version=lambda wc: str(config.get("Genome_Version", "")),
        gene_feature=_CAND_PARAMS.get("gene_feature", "gene"),
        id_attribute=_CAND_PARAMS.get("id_attribute", "ID"),
        flank_bp=_CAND_PARAMS.get("flank_bp", 0),
    log:
        "logs/11.candidate_genes/{trait}.log",
    benchmark:
        "benchmarks/11.candidate_genes/{trait}.txt",
    threads: 1,
    resources:
        **rule_resource(config, "low_resource", skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/gwasflow-python.yaml"),
    shell:
        r"""
        set -euo pipefail
        TMP_OUT="11.candidates/{wildcards.trait}/.tmp.cand.$$"
        rm -rf "$TMP_OUT"; mkdir -p "$(dirname "$TMP_OUT")"
        FT_ARG=""
        if [ -n "{input.function_table:q}" ]; then FT_ARG="--function-table {input.function_table:q}"; fi
        python3 {params.script:q} \
            --loci {input.loci:q} \
            --gff {input.gff3:q} \
            --genome-version {params.genome_version:q} \
            --gene-feature {params.gene_feature:q} \
            --id-attribute {params.id_attribute:q} \
            --flank-bp {params.flank_bp} \
            $FT_ARG \
            --out "$TMP_OUT/candidates.tsv" > {log:q} 2>&1
        mkdir -p "$(dirname {output.candidates:q})"
        mv "$TMP_OUT/candidates.tsv" {output.candidates:q}
        mv "$TMP_OUT/candidates.tsv.meta.json" {output.meta:q}
        rmdir "$TMP_OUT"
        """
