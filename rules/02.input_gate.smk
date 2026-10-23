#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# GWASFlow - 02.input_gate.smk: Input gating and genotype normalization
#
# Contracts:
# - The full set of external input files goes into rule inputs (never params, spec §5.1.1-5)
# - Gate artifact input_gate.json is a hard gate: downstream rules attach it explicitly as input
# - Rule logs stay separate from tool-owned --log-file logs
# - Temp directory + atomic move on success (design §2.1-3-4)

import os


def _genotype_input_files(cfg):
    """Return the full set of genotype input files (index for VCF/BCF; bed/bim/fam trio for BED)."""
    fmt = cfg.get("genotype_format", "vcf")
    files = []
    if fmt in ("vcf", "bcf"):
        main = cfg.get("genotype_input", "")
        files.append(main)
        index = None
        for suffix in (".csi", ".tbi"):
            candidate = main + suffix
            if os.path.exists(candidate):
                index = candidate
                break
        if index is None:
            index = main + ".csi"   # still declared when absent (gate reports the missing file explicitly)
        files.append(index)
    else:
        bfile = cfg.get("genotype_bfile", "")
        files.extend([bfile + ".bed", bfile + ".bim", bfile + ".fam"])
        if cfg.get("input_variant_map"):
            files.append(cfg["input_variant_map"])
    return files


rule verify_inputs:
    """Input gate: per-file role/component/readability checks + 3-state MD5 + table contracts + on-demand reference validation.

    Input: all external files (genotype trio/VCF+index, variant_map, samples, phenotypes, traits, covariates, expected checksums), input/reference manifests, config snapshot.
    Process: verify each file's role, components and readability; MD5 in three states (verified / mismatch / baseline_only without a baseline); validate table column contracts and reference resources required by the enabled modules only.
    Output: 01.gate/input_gate.json + md5_check.tsv + reference_check.json, and 13.summary/effective_config.yaml (Step 6 snapshot flushed here).
    """
    input:
        genotype=_genotype_input_files(config),
        samples=config["sample_csv"],
        phenotypes=config["phenotypes_tsv"],
        traits=config["traits_tsv"],
        covariates=(
            config["covariates_tsv"]
            if any(
                (t.get("numeric_covariates") or t.get("categorical_covariates"))
                for t in TRAITS.values()
            )
            else []
        ),
        expected_checksums=(
            config["expected_checksums_tsv"]
            if config.get("expected_checksums_tsv")
            else []
        ),
        input_manifest="13.summary/.input_manifest.source.json",
        reference_manifest="13.summary/.reference_manifest.source.json",
        config_snapshot="13.summary/.effective_config.source.yaml",
    output:
        gate="01.gate/input_gate.json",
        md5="01.gate/md5_check.tsv",
        reference_check="01.gate/reference_check.json",
        effective_config="13.summary/effective_config.yaml",
    params:
        script=workflow.source_path("../scripts/verify_inputs.py"),
        checksum_mode=config.get("checksum_mode", "auto"),
        threads=config["parameter"]["threads"]["verify_inputs"],
        required_covariates=(
            "true"
            if any(
                (t.get("numeric_covariates") or t.get("categorical_covariates"))
                for t in TRAITS.values()
            )
            else "false"
        ),
    log:
        "logs/02.input_gate/verify.log",
    benchmark:
        "benchmarks/02.input_gate/verify.txt",
    threads: lambda wc: config["parameter"]["threads"]["verify_inputs"],
    resources:
        **rule_resource(config, "low_resource", skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/gwasflow-python.yaml"),
    shell:
        r"""
        set -euo pipefail
        TMP_OUT="01.gate/.tmp.verify.$$"
        rm -rf "$TMP_OUT"; mkdir -p "$(dirname "$TMP_OUT")"
        COV_ARG=""
        if [ -n "{input.covariates:q}" ]; then COV_ARG="--covariates {input.covariates:q}"; fi
        CKS_ARG=""
        if [ -n "{input.expected_checksums:q}" ]; then CKS_ARG="--expected-checksums {input.expected_checksums:q}"; fi
        python3 {params.script:q} \
            --input-manifest {input.input_manifest:q} \
            --samples {input.samples:q} \
            --phenotypes {input.phenotypes:q} \
            --traits {input.traits:q} \
            $COV_ARG \
            --reference-manifest {input.reference_manifest:q} \
            $CKS_ARG \
            --checksum-mode {params.checksum_mode:q} \
            --effective-config {input.config_snapshot:q} \
            --threads {params.threads} \
            --required-covariates {params.required_covariates:q} \
            --out-dir "$TMP_OUT" > {log:q} 2>&1
        mkdir -p "$(dirname {output.gate:q})"
        mv "$TMP_OUT/input_gate.json" {output.gate:q}
        mv "$TMP_OUT/md5_check.tsv" {output.md5:q}
        mv "$TMP_OUT/reference_check.json" {output.reference_check:q}
        mkdir -p "$(dirname {output.effective_config:q})"
        mv "$TMP_OUT/effective_config.yaml" {output.effective_config:q}
        rmdir "$TMP_OUT"
        """


# VCF-route reference resources (resolved per Genome_Version)
_REFERENCE = config.get("GWAS_reference", {}).get(str(config.get("Genome_Version", "")).strip(), {})
_QC_PARAMS = config["parameter"]["qc"]
_NORM_PARAMS = config["parameter"]["normalization"]


rule normalize_genotype:
    """Genotype normalization: VCF/BCF command chain or BED cross-check into the engine format.

    Input: raw genotype files (VCF/BCF+index, or BED trio + variant_map), gate artifact, samples.csv, optional reference FASTA/FAI/chromosome map.
    Process: the VCF/BCF route runs the bcftools filter/norm/set-id/rename-chrs chain and converts to BED via plink2 (--double-id); the BED route cross-checks BIM against the provided variant_map (identity, coordinates, allele sets) and copies the trio unchanged; FAM identities are reconciled against samples.csv.
    Output: 02.genotype/raw_engine.{bed,bim,fam} + variant_map.tsv + normalize_audit.json.
    """
    input:
        genotype=_genotype_input_files(config),
        gate="01.gate/input_gate.json",
        samples=config["sample_csv"],
        genome_fasta=(_REFERENCE.get("genome_fasta") or []),
        genome_fai=(_REFERENCE.get("genome_fai") or []),
        chromosome_map=(_REFERENCE.get("chromosome_rename_map") or []),
    output:
        bed="02.genotype/raw_engine.bed",
        bim="02.genotype/raw_engine.bim",
        fam="02.genotype/raw_engine.fam",
        variant_map="02.genotype/variant_map.tsv",
        audit="02.genotype/normalize_audit.json",
    params:
        script=workflow.source_path("../scripts/normalize_genotype.py"),
        genotype_format=config.get("genotype_format", "vcf"),
        genotype_input=config.get("genotype_input", ""),
        genotype_bfile=config.get("genotype_bfile", ""),
        input_variant_map=config.get("input_variant_map", ""),
        filter_policy=_NORM_PARAMS.get("filter_policy", "pass_only"),
        dp_min=_QC_PARAMS.get("dp_min"),
        gq_min=_QC_PARAMS.get("gq_min"),
        chr_count=_QC_PARAMS.get("chr_count", 9),
        threads=config["parameter"]["threads"]["plink2"],
        memory_mb=_QC_PARAMS.get("memory_mb", 8000),
    log:
        "logs/02.input_gate/normalize.log",
    benchmark:
        "benchmarks/02.input_gate/normalize.txt",
    threads: lambda wc: config["parameter"]["threads"]["plink2"],
    resources:
        **rule_resource(config, "medium_resource", skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/gwasflow-python.yaml"),
    shell:
        r"""
        set -euo pipefail
        TMP_OUT="02.genotype/.tmp.normalize.$$"
        rm -rf "$TMP_OUT"; mkdir -p "$(dirname "$TMP_OUT")"
        REF_ARGS=""
        if [ -n "{input.genome_fasta:q}" ]; then
            REF_ARGS="--genome-fasta {input.genome_fasta:q} --fai {input.genome_fai:q} --chromosome-map {input.chromosome_map:q}"
        fi
        INPUT_ARGS="--input {params.genotype_input:q}"
        if [ -n "{params.genotype_bfile:q}" ]; then
            INPUT_ARGS="--bfile {params.genotype_bfile:q} --input-variant-map {params.input_variant_map:q}"
        fi
        DPGQ_ARGS=""
        if [ "{params.dp_min}" != "None" ] && [ -n "{params.dp_min}" ]; then
            DPGQ_ARGS="--dp-min {params.dp_min}"
        fi
        if [ "{params.gq_min}" != "None" ] && [ -n "{params.gq_min}" ]; then
            DPGQ_ARGS="$DPGQ_ARGS --gq-min {params.gq_min}"
        fi
        python3 {params.script:q} \
            --genotype-format {params.genotype_format:q} \
            $INPUT_ARGS \
            --gate {input.gate:q} \
            --samples {input.samples:q} \
            $REF_ARGS \
            --filter-policy {params.filter_policy:q} \
            $DPGQ_ARGS \
            --chr-count {params.chr_count} \
            --threads {params.threads} \
            --memory-mb {params.memory_mb} \
            --out-dir "$TMP_OUT" > {log:q} 2>&1
        mkdir -p "$(dirname {output.bed:q})"
        mv "$TMP_OUT/raw_engine.bed" {output.bed:q}
        mv "$TMP_OUT/raw_engine.bim" {output.bim:q}
        mv "$TMP_OUT/raw_engine.fam" {output.fam:q}
        mv "$TMP_OUT/variant_map.tsv" {output.variant_map:q}
        mv "$TMP_OUT/normalize_audit.json" {output.audit:q}
        rm -rf "$TMP_OUT"
        """
