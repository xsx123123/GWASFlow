#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GWASFlow → BioReport 报告配置导出（M6 衔接）。

从 GWASFlow 工作区产物自动生成 bioreport 的 report_config.json：
- project 段：13.summary/project_summary.json（project_name/species/genome_version/pipeline_version）
- traits：08.results/{trait}/metadata.json 逐性状统计（analysis.gwas.traits）
- 软件版本：--software-versions JSON 注入 > GEMMA help 抽取 + provenance 兜底
- input_files：QC 审计/关联/loci/PCA/亲缘/染色体表（相对 workflow 根，path_base=workflow 根）

错误码：0=成功；2=输入/结构错误。
"""

import argparse
import json
import os
import re
import sys
from pathlib import Path

EXIT_OK = 0
EXIT_INPUT_ERROR = 2


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="export_report_config.py",
        description="GWASFlow 工作区 → BioReport report_config.json（报告配置自动导出）",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--workflow-dir", required=True, help="GWASFlow 工作区根（01.gate..14.data_deliver 所在目录）")
    parser.add_argument("--out", required=True, help="report_config.json 输出路径（建议 15.report/ 下）")
    parser.add_argument("--project-id", default="", help="报告项目编号（默认 PRJ-<project_name>）")
    parser.add_argument("--client", default="", help="委托单位/客户名（报告封面展示）")
    parser.add_argument("--chromosome-table", default="",
                        help="染色体表覆盖（默认读 effective_config 的 GWAS_reference）")
    parser.add_argument("--software-versions", default="",
                        help='软件版本 JSON 数组（如 \'[{"Software":"GEMMA","Version":"0.98.5"}]\'）')
    parser.add_argument("--output-format", default="html,pdf", help="报告输出格式（逗号分隔：html,pdf）")
    return parser.parse_args(argv)


def fail(message, code=EXIT_INPUT_ERROR):
    print(f"[export_report_config ERROR] {message}", file=sys.stderr)
    sys.exit(code)


def main(argv=None):
    args = parse_args(argv)
    workflow = Path(args.workflow_dir).resolve()
    if not workflow.is_dir():
        fail(f"工作区目录不存在: {workflow}")

    summary_path = workflow / "13.summary" / "project_summary.json"
    if not summary_path.exists():
        fail(f"缺少项目摘要（13.summary/project_summary.json）: {summary_path}")
    with open(summary_path, encoding="utf-8") as handle:
        summary = json.load(handle)

    project = summary.get("project", {})
    counts = summary.get("counts", {})

    # ---- 性状清单（08.results/{trait}/metadata.json）----
    traits = []
    results_dir = workflow / "08.results"
    if results_dir.is_dir():
        for child in sorted(results_dir.iterdir()):
            if (child / "metadata.json").is_file():
                traits.append(child.name)

    # ---- 染色体表：参数覆盖 > effective_config ----
    chromosome_table = args.chromosome_table
    if not chromosome_table:
        effective_path = workflow / "13.summary" / "effective_config.yaml"
        if effective_path.exists():
            try:
                import yaml
                with open(effective_path, encoding="utf-8") as handle:
                    effective = yaml.safe_load(handle) or {}
                version = str(effective.get("Genome_Version", ""))
                chromosome_table = str(
                    (effective.get("GWAS_reference", {}) or {}).get(version, {}).get("chromosome_table", ""))
            except Exception as error:  # 染色体表缺失时报告曼哈顿图无法构建基因组坐标
                fail(f"effective_config.yaml 解析失败: {error}")
    if not chromosome_table:
        fail("染色体表未解析（--chromosome-table 或 GWAS_reference.<version>.chromosome_table）")

    # ---- 亲缘矩阵名（centered=cXX / standardized=sXX）----
    kinship_matrix = "06.kinship/kinship.cXX.txt"
    effective_path = workflow / "13.summary" / "effective_config.yaml"
    if effective_path.exists():
        try:
            import yaml
            with open(effective_path, encoding="utf-8") as handle:
                effective = yaml.safe_load(handle) or {}
            k_type = ((effective.get("parameter", {}) or {}).get("kinship", {}) or {}).get("type", "centered")
            kinship_matrix = "06.kinship/kinship.cXX.txt" if k_type == "centered" else "06.kinship/kinship.sXX.txt"
        except Exception:
            pass

    # ---- 软件版本：注入 > GEMMA help 抽取 + provenance 兜底 ----
    if args.software_versions:
        try:
            software_versions = json.loads(args.software_versions)
        except json.JSONDecodeError as error:
            fail(f"--software-versions JSON 解析失败: {error}")
    else:
        software_versions = []
        gemma_help = None
        for help_path in sorted((workflow / "07.association").glob("*/gemma.help.txt")):
            gemma_help = help_path
            break
        if gemma_help and gemma_help.exists():
            match = re.search(r"GEMMA\s+([0-9.]+)", gemma_help.read_text(encoding="utf-8", errors="ignore"))
            if match:
                software_versions.append({"Software": "GEMMA", "Version": match.group(1)})
        provenance = summary.get("provenance", {}) or {}
        for name, version in (provenance.get("software_versions") or {}).items():
            if name == "pipeline" or "占位" in str(version):
                continue
            software_versions.append({"Software": name, "Version": str(version)})
        software_versions.append({"Software": "GWASFlow", "Version": str(project.get("pipeline_version", ""))})

    # ---- input_files（相对 workflow 根；bioreport 以 path_base 解析）----
    def rel(path):
        return str(Path(path).relative_to(workflow))

    input_files = {
        "project_summary_file": rel(summary_path),
        "chromosome_table_file": chromosome_table if os.path.isabs(chromosome_table)
        else chromosome_table,
        "gwas_dir": "08.results",
        "loci_dir": "09.loci",
    }
    for key, candidate in (
        ("qc_final_stats_file", workflow / "03.genotype_qc" / "final_stats.json"),
        ("qc_audit_file", workflow / "03.genotype_qc" / "qc_audit.json"),
        ("pca_eigenvec_file", workflow / "04.structure" / "pca.eigenvec"),
        ("kinship_matrix_file", workflow / kinship_matrix),
    ):
        if candidate.exists():
            input_files[key] = rel(candidate)

    config = {
        "schema_version": "2.0.0",
        "path_base": str(workflow),
        "input_files": input_files,
        "software_versions": software_versions,
        "project": {
            "project_id": args.project_id or f"PRJ-{project.get('project_name', 'GWASFLOW')}",
            "client": args.client,
            "species": project.get("species", ""),
            "genome_version": project.get("genome_version", ""),
            "pipeline_version": project.get("pipeline_version", ""),
        },
        "analysis": {
            "gwas": {
                "traits": traits,
                "manhattan": {"point_budget": 200000},
                "significant": {"top_n": 20},
                "qq": {"max_points": 50000},
                "kinship": {"heatmap_max_samples": 194},
            }
        },
        "stats": {
            "total_samples": counts.get("n_samples_final"),
            "trait_count": len(traits),
        },
        "ai": {"enabled": False},
        "report": {
            "type": "gwas",
            "layout": "single",
            "output_format": [f.strip() for f in args.output_format.split(",") if f.strip()],
            "modules": ["gwas.index", "gwas.workflow", "gwas.qc", "gwas.results"],
            "title": "GWAS 全基因组关联分析报告",
            "title_en": "GWAS Analysis Report",
        },
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as handle:
        json.dump(config, handle, ensure_ascii=False, indent=2)

    missing = [key for key in ("qc_final_stats_file", "qc_audit_file") if key not in input_files]
    if missing:
        print(f"[export_report_config WARNING] 以下输入缺失（对应报告章节将被跳过）: {missing}", file=sys.stderr)
    print(f"[export_report_config] {out}: {len(traits)} 性状 {traits}，"
          f"软件版本 {len(software_versions)} 项，染色体表 {chromosome_table}")
    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
