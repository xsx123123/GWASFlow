#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# GWASFlow - 06.kinship.smk: Cohort-wide kinship matrix (computation + validation)
#
# Cohort-wide K (the all-ones phenotype forces full coverage); the sample set of K is always the full FAM set.
# Dense K memory scales as N^2: gwas_kinship_large high-memory profile (no single fixed mem may promise to fit all cohorts).
# Rule logs stay separate from GEMMA's own logs.

_KINSHIP_PARAMS = config["parameter"]["kinship"]
_K_MATRIX = "kinship.cXX.txt" if (_KINSHIP_PARAMS.get("type", "centered") == "centered") else "kinship.sXX.txt"


def _first_trait_pheno_all():
    """Any trait's all-ones phenotype file (K is shared cohort-wide)."""
    if TRAITS:
        first_trait = sorted(TRAITS.keys())[0]
        return f"05.phenotype/{first_trait}/kinship_all.pheno.txt"
    return []


rule compute_kinship:
    """Cohort-wide kinship matrix via GEMMA -gk 1|2 (-p all-ones phenotype, -bfile cohort).

    Input: 03.genotype_qc/cohort.{bed,bim,fam} + any trait's kinship_all.pheno.txt + gate artifact.
    Process: the wrapper script builds and runs GEMMA -gk 1 (centered) or -gk 2 (standardized) with the all-ones phenotype so every cohort sample enters K; the GEMMA command and run info are recorded for reconciliation.
    Output: 06.kinship/kinship.cXX.txt (or .sXX.txt) + kinship.log.txt + sample_order.tsv.
    """
    input:
        bed="03.genotype_qc/cohort.bed",
        bim="03.genotype_qc/cohort.bim",
        fam="03.genotype_qc/cohort.fam",
        pheno_all=_first_trait_pheno_all(),
        gate="01.gate/input_gate.json",
    output:
        matrix=f"06.kinship/{_K_MATRIX}",
        gemma_log="06.kinship/kinship.log.txt",
        sample_order="06.kinship/sample_order.tsv",
    params:
        script=workflow.source_path("../scripts/compute_kinship.py"),
        kinship_type=_KINSHIP_PARAMS.get("type", "centered"),
        threads=config["parameter"]["threads"]["gemma"],
    log:
        "logs/06.kinship/compute_kinship.log",
    benchmark:
        "benchmarks/06.kinship/compute_kinship.txt",
    threads: lambda wc: config["parameter"]["threads"]["gemma"],
    resources:
        **rule_resource(config, "gwas_kinship_large", skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/gemma.yaml"),
    shell:
        r"""
        set -euo pipefail
        TMP_OUT="06.kinship/.tmp.kinship.$$"
        rm -rf "$TMP_OUT"; mkdir -p "$(dirname "$TMP_OUT")"
        python3 {params.script:q} \
            --bfile 03.genotype_qc/cohort \
            --pheno-all-ones {input.pheno_all:q} \
            --fam {input.fam:q} \
            --mode genome_wide \
            --type {params.kinship_type:q} \
            --threads {params.threads} \
            --out-dir "$TMP_OUT" > {log:q} 2>&1
        mkdir -p "$(dirname {output.matrix:q})"
        mv "$TMP_OUT/{_K_MATRIX}" {output.matrix:q}
        mv "$TMP_OUT/kinship.log.txt" {output.gemma_log:q}
        mv "$TMP_OUT/sample_order.tsv" {output.sample_order:q}
        cp "$TMP_OUT/gemma.run.log" "06.kinship/gemma.run.log" || true
        cp "$TMP_OUT/kinship_run_info.json" "06.kinship/kinship_run_info.json" || true
        rm -rf "$TMP_OUT"
        """


rule validate_kinship:
    """Kinship matrix validation: dimensions / symmetry / row-column identity / finiteness / PSD (tolerating numerical error only).

    Input: 06.kinship kinship matrix + sample_order.tsv + cohort.fam.
    Process: check matrix dimensions against the FAM, symmetry within tolerance, row/column identity agreement, finiteness of all entries, and positive semi-definiteness within the eigenvalue tolerance.
    Output: 06.kinship/kinship_validation.json.
    """
    input:
        matrix=f"06.kinship/{_K_MATRIX}",
        matrix_ids="06.kinship/sample_order.tsv",
        fam="03.genotype_qc/cohort.fam",
    output:
        validation="06.kinship/kinship_validation.json",
    params:
        script=workflow.source_path("../scripts/validate_kinship.py"),
        symmetry_tolerance=_KINSHIP_PARAMS.get("symmetry_tolerance", 1e-8),
        eigenvalue_tolerance=_KINSHIP_PARAMS.get("eigenvalue_tolerance", 1e-8),
    log:
        "logs/06.kinship/validate_kinship.log",
    benchmark:
        "benchmarks/06.kinship/validate_kinship.txt",
    threads: 1,
    resources:
        **rule_resource(config, "low_resource", skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/gemma.yaml"),
    shell:
        r"""
        set -euo pipefail
        python3 {params.script:q} \
            --matrix {input.matrix:q} \
            --matrix-ids {input.matrix_ids:q} \
            --fam {input.fam:q} \
            --symmetry-tolerance {params.symmetry_tolerance} \
            --eigenvalue-tolerance {params.eigenvalue_tolerance} \
            --out {output.validation:q} > {log:q} 2>&1
        """
