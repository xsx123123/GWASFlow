#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# GWASFlow - 15.Report.smk: 条件报告（M6 边界）
#
# report 默认 false；report=true 且渲染器未验收（M6）→ run 段显式报错退出 2，
# 指明 GWAS 渲染器属后续工作包/条件接入；BioReport 不具备时交付结构化摘要+静态图表。
# deliver/report 两开关独立、不互相隐式关闭。

rule render_report:
    """GWAS 报告渲染（条件实现：渲染器验收前显式报错）。"""
    input:
        project_summary="13.summary/project_summary.json",
        plot_manifest="12.plots/plot_manifest.json",
        manifest="14.data_deliver/delivery_manifest.json",
    output:
        report_manifest="15.report/report_manifest.json",
    log:
        "logs/15.Report/render_report.log",
    benchmark:
        "benchmarks/15.Report/render_report.txt",
    threads: 1,
    resources:
        **rule_resource(config, "low_resource", skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/gwasflow-python.yaml"),
    shell:
        r"""
        set -euo pipefail
        echo "[render_report ERROR] GWAS 报告渲染器属后续工作包/条件接入（M6：BioReport 本地不可用，" \
             "渲染器未验收）。当前能力边界：交付结构化摘要（project_summary.json）" \
             "与静态图表（12.plots/）。如需启用报告，请先完成 BioReport GWAS 渲染器" \
             "验收（独立工作包）。" > {log:q}
        cat {log:q} >&2
        exit 2
        """
