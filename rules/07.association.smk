#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# GWASFlow - 07.association.smk: Per-trait MAC + GEMMA LMM (per trait, D5)
#
# No cross-trait dependency edges ({trait} wildcard expands independently, spec §4.1-1).
# GEMMA uses the full cohort bfile (D10: comparable with independent reference commands at identical parameters).
# Rule logs stay separate from GEMMA's own --log-file.

_ASSOC_PARAMS = config["parameter"]["association"]
_KINSHIP_PARAMS = config["parameter"]["kinship"]
_K_MATRIX = "kinship.cXX.txt" if (_KINSHIP_PARAMS.get("type", "centered") == "centered") else "kinship.sXX.txt"


def _assoc_output_files(test):
    mapping = {
        "wald": ["association.assoc.txt"],
        "lrt": ["association.lrt.txt"],
        "score": ["association.score.txt"],
        "all": ["association.assoc.txt", "association.lrt.txt", "association.score.txt"],
    }
    return mapping[test]


_RUN_GEMMA_OUTPUTS = {
    "gemma_log": "07.association/{trait}/association.log.txt",
    "help_text": "07.association/{trait}/gemma.help.txt",
    "run_meta": "07.association/{trait}/run_metadata.json",
}
for _f in _assoc_output_files(_ASSOC_PARAMS.get("primary_test", "wald")):
    _RUN_GEMMA_OUTPUTS[_f.replace("association.", "").replace(".txt", "")] = (
        f"07.association/{{trait}}/{_f}"
    )


rule trait_mac:
    """Exact MAC on the trait's effective samples (plink2 --keep --freq, D4).

    Input: 03.genotype_qc/cohort.{bed,bim,fam} + 05.phenotype/{trait}/keep.tsv + gate artifact.
    Process: run plink2 --keep <keep.tsv> --freq on the cohort to obtain exact allele counts, then compute per-variant minor allele counts restricted to the trait's effective samples.
    Output: 07.association/{trait}/trait_mac.tsv + trait_mac_summary.json.
    """
    input:
        bed="03.genotype_qc/cohort.bed",
        bim="03.genotype_qc/cohort.bim",
        fam="03.genotype_qc/cohort.fam",
        keep="05.phenotype/{trait}/keep.tsv",
        gate="01.gate/input_gate.json",
    output:
        mac="07.association/{trait}/trait_mac.tsv",
        summary="07.association/{trait}/trait_mac_summary.json",
    params:
        script=workflow.source_path("../scripts/compute_trait_mac.py"),
        mac_min=config["parameter"]["qc"].get("mac", 5),
        threads=config["parameter"]["threads"]["plink2"],
    log:
        "logs/07.association/{trait}_mac.log",
    benchmark:
        "benchmarks/07.association/{trait}_mac.txt",
    threads: lambda wc: config["parameter"]["threads"]["plink2"],
    resources:
        **rule_resource(config, "medium_resource", skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/plink2.yaml"),
    shell:
        r"""
        set -euo pipefail
        TMP_OUT="07.association/{wildcards.trait}/.tmp.mac.$$"
        rm -rf "$TMP_OUT"; mkdir -p "$(dirname "$TMP_OUT")"
        python3 {params.script:q} \
            --bfile 03.genotype_qc/cohort \
            --fam {input.fam:q} \
            --keep {input.keep:q} \
            --mac-min {params.mac_min} \
            --threads {params.threads} \
            --out-dir "$TMP_OUT" > {log:q} 2>&1
        mkdir -p "$(dirname {output.mac:q})"
        mv "$TMP_OUT/trait_mac.tsv" {output.mac:q}
        mv "$TMP_OUT/trait_mac_summary.json" {output.summary:q}
        rm -rf "$TMP_OUT"
        """


rule run_gemma:
    """GEMMA LMM association test (per trait; command-array construction + run_metadata five-way binding).

    Input: cohort trio + trait GEMMA phenotype/covariates + trait metadata/sample order + kinship matrix/order/validation + trait MAC + gate artifact.
    Process: the wrapper script builds the GEMMA command (-lmm <test>, -c covariates, -k K, explicit -maf/-miss/-hwe/-r2) and binds five artifact classes (inputs, trait, kinship, assoc, binary) into run_metadata.json; --keep-going lets GEMMA finish despite per-variant failures.
    Output: 07.association/{trait}/association.assoc.txt (+ .lrt/.score when test=all) + association.log.txt + gemma.help.txt + run_metadata.json.
    """
    input:
        bed="03.genotype_qc/cohort.bed",
        bim="03.genotype_qc/cohort.bim",
        fam="03.genotype_qc/cohort.fam",
        pheno="05.phenotype/{trait}/phenotype.gemma.txt",
        covar="05.phenotype/{trait}/covariates.gemma.txt",
        trait_metadata="05.phenotype/{trait}/metadata.json",
        trait_order="05.phenotype/{trait}/sample_order.tsv",
        kinship=f"06.kinship/{_K_MATRIX}",
        kinship_order="06.kinship/sample_order.tsv",
        kinship_validation="06.kinship/kinship_validation.json",
        trait_mac="07.association/{trait}/trait_mac.tsv",
        gate="01.gate/input_gate.json",
    output:
        **_RUN_GEMMA_OUTPUTS,
    params:
        script=workflow.source_path("../scripts/run_gemma.py"),
        trait=lambda wc: wc.trait,
        trait_dir=lambda wc: f"05.phenotype/{wc.trait}",
        kinship_type=_KINSHIP_PARAMS.get("type", "centered"),
        test=_ASSOC_PARAMS.get("primary_test", "wald"),
        maf=_ASSOC_PARAMS.get("maf", 0.05),
        miss=_ASSOC_PARAMS.get("miss", 0.10),
        threads=config["parameter"]["threads"]["gemma"],
        resource_profile="high_resource",
    log:
        "logs/07.association/{trait}.log",
    benchmark:
        "benchmarks/07.association/{trait}.txt",
    threads: lambda wc: config["parameter"]["threads"]["gemma"],
    resources:
        **rule_resource(config, "high_resource", skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/gemma.yaml"),
    shell:
        r"""
        set -euo pipefail
        TMP_OUT="07.association/{wildcards.trait}/.tmp.gemma.$$"
        rm -rf "$TMP_OUT"; mkdir -p "$(dirname "$TMP_OUT")"
        python3 {params.script:q} \
            --bfile 03.genotype_qc/cohort \
            --trait-dir {params.trait_dir:q} \
            --trait {params.trait:q} \
            --kinship {input.kinship:q} \
            --kinship-sample-order {input.kinship_order:q} \
            --kinship-type {params.kinship_type:q} \
            --trait-mac {input.trait_mac:q} \
            --test {params.test:q} \
            --maf {params.maf} \
            --miss {params.miss} \
            --threads {params.threads} \
            --resource-profile {params.resource_profile:q} \
            --keep-going \
            --out-dir "$TMP_OUT" > {log:q} 2>&1
        mkdir -p "$(dirname {output.run_meta:q})"
        mv "$TMP_OUT/association.log.txt" {output.gemma_log:q}
        mv "$TMP_OUT/gemma.help.txt" {output.help_text:q}
        mv "$TMP_OUT/run_metadata.json" {output.run_meta:q}
        for f in association.assoc.txt association.lrt.txt association.score.txt; do
            if [ -f "$TMP_OUT/$f" ]; then
                mv "$TMP_OUT/$f" "07.association/{wildcards.trait}/$f"
            fi
        done
        rm -rf "$TMP_OUT"
        """
