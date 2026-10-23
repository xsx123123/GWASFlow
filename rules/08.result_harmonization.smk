#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# GWASFlow - 08.result_harmonization.smk: Harmonized association statistics (five artifacts per trait)
#
# variant_map is supplied directly from NORM (not forwarded through QC): the reconciliation baseline is the original identity at normalization time.

_ASSOC_PARAMS = config["parameter"]["association"]
_MT_PARAMS = config["parameter"]["multiple_testing"]
_TRAIT_MODEL = "gemma_lmm"


def _assoc_input(wildcards):
    test = _ASSOC_PARAMS.get("primary_test", "wald")
    mapping = {
        "wald": "07.association/{trait}/association.assoc.txt",
        "lrt": "07.association/{trait}/association.lrt.txt",
        "score": "07.association/{trait}/association.score.txt",
        "all": "07.association/{trait}/association.assoc.txt",
    }
    return mapping[test].format(trait=wildcards.trait)


rule harmonize:
    """Harmonized association statistics: strict identity validation + direction conversion + tiny-P handling + Bonferroni + exact MAC.

    Input: raw GEMMA assoc output + 02.genotype/variant_map.tsv + trait sample metadata + run_metadata + trait MAC + gate artifact.
    Process: strong identity validation (rs -> variant_map lookup, chr/pos match, allele-set match; structural errors abort, no silent skips); effect-direction conversion per target_effect_allele (reported keeps GEMMA allele1, alt flips beta/eaf); Decimal-precision tiny P (p=0 recorded as-is, no fabricated values); Bonferroni threshold = alpha/(family_size x M) with M = valid P count; exact MAC attached per variant; unreported variants reconciled against variant_map.
    Output: 08.results/{trait}/gwas.tsv + significant.tsv + invalid_variants.tsv + unreported_variants.tsv + metadata.json.
    """
    input:
        assoc=lambda wc: _assoc_input(wc),
        variant_map="02.genotype/variant_map.tsv",
        sample_metadata="05.phenotype/{trait}/metadata.json",
        run_metadata="07.association/{trait}/run_metadata.json",
        trait_mac="07.association/{trait}/trait_mac.tsv",
        gate="01.gate/input_gate.json",
    output:
        gwas="08.results/{trait}/gwas.tsv",
        significant="08.results/{trait}/significant.tsv",
        invalid="08.results/{trait}/invalid_variants.tsv",
        unreported="08.results/{trait}/unreported_variants.tsv",
        metadata="08.results/{trait}/metadata.json",
    params:
        script=workflow.source_path("../scripts/harmonize_gemma.py"),
        trait=lambda wc: wc.trait,
        analysis_id=lambda wc: f"{wc.trait}_gemma_primary_001",
        genome_version=lambda wc: str(config.get("Genome_Version", "")),
        effect_scale=lambda wc: TRAITS[wc.trait]["effect_scale"],
        p_column=lambda wc: {
            "wald": "p_wald", "lrt": "p_lrt", "score": "p_score",
            "all": "p_wald",
        }[_ASSOC_PARAMS.get("primary_test", "wald")],
        target_effect_allele=config.get("target_effect_allele", "reported"),
        alpha=str(_MT_PARAMS.get("alpha", "0.05")),
        family_size=_MT_PARAMS.get("family_size", 1),
    log:
        "logs/08.result_harmonization/{trait}.log",
    benchmark:
        "benchmarks/08.result_harmonization/{trait}.txt",
    threads: 1,
    resources:
        **rule_resource(config, "low_resource", skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/gwasflow-python.yaml"),
    shell:
        r"""
        set -euo pipefail
        TMP_OUT="08.results/{wildcards.trait}/.tmp.harm.$$"
        rm -rf "$TMP_OUT"; mkdir -p "$(dirname "$TMP_OUT")"
        set +e
        python3 {params.script:q} \
            --assoc {input.assoc:q} \
            --variant-map {input.variant_map:q} \
            --sample-metadata {input.sample_metadata:q} \
            --run-metadata {input.run_metadata:q} \
            --trait-mac {input.trait_mac:q} \
            --analysis-id {params.analysis_id:q} \
            --genome-version {params.genome_version:q} \
            --effect-scale {params.effect_scale:q} \
            --trait {params.trait:q} \
            --p-column {params.p_column:q} \
            --target-effect-allele {params.target_effect_allele:q} \
            --alpha {params.alpha:q} \
            --family-size {params.family_size} \
            --out-dir "$TMP_OUT" > {log:q} 2>&1
        HARM_EXIT=$?
        set -e
        if [ "$HARM_EXIT" -ne 0 ] && [ "$HARM_EXIT" -ne 3 ]; then
            exit "$HARM_EXIT"
        fi
        mkdir -p "$(dirname {output.gwas:q})"
        mv "$TMP_OUT/gwas.tsv" {output.gwas:q}
        mv "$TMP_OUT/significant.tsv" {output.significant:q}
        mv "$TMP_OUT/invalid_variants.tsv" {output.invalid:q}
        mv "$TMP_OUT/unreported_variants.tsv" {output.unreported:q}
        mv "$TMP_OUT/metadata.json" {output.metadata:q}
        rm -rf "$TMP_OUT"
        exit "$HARM_EXIT"
        """
