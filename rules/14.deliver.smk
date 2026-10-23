#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# GWASFlow - 14.deliver.smk: Copy delivery (manifest trio)
#
# input reuses the snakefile-level ALL_TARGETS (single-call-site constraint, framework spec §5.6-2);
# delivery is based on the explicit artifact_index set (recursive globbing forbidden); the default copy mode stays independently movable.

rule deliver_results:
    """Copy delivery by the explicit artifact_index + per-file MD5 + manifest trio.

    Input: artifact_index.tsv + project_summary.json (+ all analysis/summary targets via ALL_TARGETS).
    Process: copy each artifact listed in the explicit index into the delivery tree; compute per-file MD5; write the delivery manifest trio; mirror the delivery directory to the external deliver root when configured.
    Output: 14.data_deliver/delivery_manifest.json + delivery_manifest.md5 + delivery_details.log (+ analysis_results/ tree).
    """
    input:
        analysis=lambda wc: list(ALL_TARGETS.plan_analysis),
        summary=lambda wc: list(ALL_TARGETS.plan_summary),
        artifact_index="13.summary/artifact_index.tsv",
        project_summary="13.summary/project_summary.json",
    output:
        manifest="14.data_deliver/delivery_manifest.json",
        manifest_md5="14.data_deliver/delivery_manifest.md5",
        details="14.data_deliver/delivery_details.log",
    params:
        script=workflow.source_path("../scripts/deliver_results.py"),
        source_dir=config["workflow"],
        target_dir="14.data_deliver/analysis_results",
        out_dir="14.data_deliver",
        deliver_root=config["data_deliver"],
        threads=config["parameter"]["threads"]["deliver"],
    log:
        "logs/14.deliver/deliver.log",
    benchmark:
        "benchmarks/14.deliver/deliver.txt",
    threads: lambda wc: config["parameter"]["threads"]["deliver"],
    resources:
        **rule_resource(config, "low_resource", skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/gwasflow-python.yaml"),
    shell:
        r"""
        set -euo pipefail
        python3 {params.script:q} \
            --artifact-index {input.artifact_index:q} \
            --source-dir {params.source_dir:q} \
            --target-dir {params.target_dir:q} \
            --mode copy \
            --threads {params.threads} \
            --out-dir {params.out_dir:q} > {log:q} 2>&1
        # Mirror to the external delivery root (optional; keeps the delivery tree independently movable)
        if [ -d "{params.out_dir:q}" ] && [ -n "{params.deliver_root:q}" ] && [ "{params.deliver_root:q}" != "{params.out_dir:q}" ]; then
            mkdir -p "{params.deliver_root:q}"
            cp -R "{params.out_dir:q}/." "{params.deliver_root:q}/" >> {log:q} 2>&1 || true
        fi
        """
