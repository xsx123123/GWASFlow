#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# GWASFlow - 15.Report.smk: BioReport rendering (M6 integration)
#
# report defaults to false. With report: true:
# - exports report_config.json from the workflow outputs (export_report_config.py),
# - renders bilingual (zh/en) HTML + PDF via the BioReport CLI when available,
# - fails loudly with exit 2 (explicit boundary, never silently skipped) when the
#   BioReport engine is not installed on PATH.
# The deliver/report switches are independent and never implicitly disable each other.

_KINSHIP_TYPE = config["parameter"]["kinship"].get("type", "centered")
_KINSHIP_MATRIX = "kinship.cXX.txt" if _KINSHIP_TYPE == "centered" else "kinship.sXX.txt"


rule render_report:
    """GWAS report rendering via BioReport (bilingual zh/en, HTML + PDF; explicit error when the engine is absent).

    Input: project_summary.json + qc_audit/final_stats + per-trait association metadata/gwas + loci tables + PCA + kinship + effective_config.yaml + plot manifest + delivery manifest.
    Process: export a BioReport report_config.json from the GWASFlow outputs (traits, counts, chromosome table, software versions auto-derived); invoke `bioreport` to render the gwas.* module set in both languages and both formats; emit a report_manifest.json listing the delivered artifacts. Missing BioReport engine exits 2 with an explicit message (M6 boundary).
    Output: 15.report/report_manifest.json (+ report_config.json and the rendered report_output/ tree).
    """
    input:
        project_summary="13.summary/project_summary.json",
        effective_config="13.summary/effective_config.yaml",
        qc_audit="03.genotype_qc/qc_audit.json",
        qc_final_stats="03.genotype_qc/final_stats.json",
        pca_eigenvec="04.structure/pca.eigenvec",
        kinship=f"06.kinship/{_KINSHIP_MATRIX}",
        assoc_metadata=lambda wc: [f"08.results/{t}/metadata.json" for t in TRAITS.keys()],
        gwas_tsv=lambda wc: [f"08.results/{t}/gwas.tsv" for t in TRAITS.keys()],
        loci=lambda wc: [f"09.loci/{t}/loci.tsv" for t in TRAITS.keys()],
        plot_manifest="12.plots/plot_manifest.json",
        delivery_manifest="14.data_deliver/delivery_manifest.json",
    output:
        report_manifest="15.report/report_manifest.json",
    params:
        export_script=workflow.source_path("../scripts/export_report_config.py"),
        manifest_script=workflow.source_path("../scripts/build_report_manifest.py"),
        workflow_dir=config["workflow"],
        output_format=config.get("report_output_format", "html,pdf"),
        project_id=config.get("project_name", "gwasflow_project"),
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
        mkdir -p 15.report
        if ! command -v bioreport >/dev/null 2>&1; then
            echo "[render_report ERROR] BioReport engine not found on PATH. " \
                 "Install bioreport (pip install -e /path/to/bioreport) to enable report rendering; " \
                 "without it the deliverable stays the structured summary plus static/interactive plots " \
                 "(M6 boundary, explicit failure, never silently skipped)." | tee {log:q} >&2
            exit 2
        fi
        python3 {params.export_script:q} \
            --workflow-dir {params.workflow_dir:q} \
            --project-id {params.project_id:q} \
            --output-format {params.output_format:q} \
            --out 15.report/report_config.json >> {log:q} 2>&1
        bioreport -c 15.report/report_config.json \
            -o 15.report/report_output >> {log:q} 2>&1
        python3 {params.manifest_script:q} \
            --report-output-dir 15.report/report_output \
            --config 15.report/report_config.json \
            --out {output.report_manifest:q} >> {log:q} 2>&1
        """
