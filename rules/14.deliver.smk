#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# GWASFlow - 14.deliver.smk: 复制交付（manifest 三件套）
#
# input 复用 snakefile 层 ALL_TARGETS（唯一调用点约束，框架规范 §5.6-2）；
# 交付基于 artifact_index 显式集合（禁止递归 glob）；默认 copy 可独立移动。

rule deliver_results:
    """按显式 artifact_index 复制交付 + 逐文件 MD5 + manifest 三件套。"""
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
        # 外部交付根目录镜像（可选；保持交付目录可独立移动语义）
        if [ -d "{params.out_dir:q}" ] && [ -n "{params.deliver_root:q}" ] && [ "{params.deliver_root:q}" != "{params.out_dir:q}" ]; then
            mkdir -p "{params.deliver_root:q}"
            cp -R "{params.out_dir:q}/." "{params.deliver_root:q}/" >> {log:q} 2>&1 || true
        fi
        """
