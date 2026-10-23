#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# GWASFlow - 13.export_summary.smk: 项目摘要与显式产物索引

import json
#
# partial 不写 complete；only_qc 不引用未运行模块；artifact_index 为显式集合。


def _export_plan_json(wildcards):
    """运行计划 JSON（module_status / analysis_targets / trait metadata 路径）。"""
    trait_metadata = {}
    prepare_metadata = {}
    if MODULE_STATUS["result_harmonization"].enabled:
        for trait in TRAITS.keys():
            trait_metadata[trait] = f"08.results/{trait}/metadata.json"
    if MODULE_STATUS["phenotype_qc"].enabled:
        for trait in TRAITS.keys():
            prepare_metadata[trait] = f"05.phenotype/{trait}/metadata.json"
    return json.dumps({
        "module_status": {
            mod: {
                "enabled": status.enabled,
                "auto_enable_reasons": list(status.auto_enable_reasons),
            }
            for mod, status in MODULE_STATUS.items()
        },
        "analysis_targets": list(ALL_TARGETS.plan_analysis),
        "trait_metadata": trait_metadata,
        "prepare_metadata": prepare_metadata,
    })


rule export_summary:
    """project_summary.json（七类审计）+ artifact_index.tsv（显式集合）。"""
    input:
        analysis=lambda wc: list(ALL_TARGETS.plan_analysis),
        effective_config="13.summary/effective_config.yaml",
        qc_audit="03.genotype_qc/qc_audit.json",
        gate="01.gate/input_gate.json",
        cohort_fam="03.genotype_qc/cohort.fam",
        cohort_bim="03.genotype_qc/cohort.bim",
        sample_order=lambda wc: (
            "06.kinship/sample_order.tsv"
            if MODULE_STATUS["kinship_analysis"].enabled
            else "03.genotype_qc/cohort.fam"
        ),
    output:
        summary="13.summary/project_summary.json",
        artifact_index="13.summary/artifact_index.tsv",
    params:
        script=workflow.source_path("../scripts/export_project_summary.py"),
        plan_b64=lambda wc: __import__("base64").b64encode(
            _export_plan_json(wc).encode()
        ).decode(),
        pipeline_version=config.get("pipeline_version", "0.1.0"),
        genome_version=lambda wc: str(config.get("Genome_Version", "")),
        species=config.get("species", "unspecified"),
        project_name=config.get("project_name", "gwasflow_project"),
        seed=config.get("seed", 20261023),
        only_qc=bool(config.get("only_qc", False)),
    log:
        "logs/13.export_summary/export.log",
    benchmark:
        "benchmarks/13.export_summary/export.txt",
    threads: 1,
    resources:
        **rule_resource(config, "low_resource", skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/gwasflow-python.yaml"),
    shell:
        r"""
        set -euo pipefail
        TMP_DIR=".tmp.export.$$"
        rm -rf "$TMP_DIR"; mkdir -p "$TMP_DIR"
        echo {params.plan_b64:q} | base64 -d > "$TMP_DIR/plan.json"
        python3 {params.script:q} \
            --effective-config {input.effective_config:q} \
            --plan-json "$TMP_DIR/plan.json" \
            --qc-audit {input.qc_audit:q} \
            --gate-json {input.gate:q} \
            --cohort-fam {input.cohort_fam:q} \
            --cohort-bim {input.cohort_bim:q} \
            --sample-order {input.sample_order:q} \
            --pipeline-version {params.pipeline_version:q} \
            --genome-version {params.genome_version:q} \
            --species {params.species:q} \
            --project-name {params.project_name:q} \
            --seed {params.seed} \
            --out {output.summary:q} \
            --artifact-index-out {output.artifact_index:q} > {log:q} 2>&1
        rm -rf "$TMP_DIR"
        """
