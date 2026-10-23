#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GWASFlow 项目摘要与产物索引（参数手册 §9.8）。

project_summary.json（schema: gwasflow.project.v1）：
- 模块状态（enabled/disabled/自动启用原因链/运行状态）
- 分析索引（每性状状态与产物路径）
- 计数（样本/标记/各阶段流失/有效检验数 M/未报告位点数）
- 参数与版本（生效参数、软件版本、参考与注释版本、代码提交、种子）
- 七类审计信息（输入摘要/配置快照/代码提交/软件环境/参考版本/样本顺序/产物摘要）
- 警告聚合与未验证项清单（M3–M6 降级记录）

红线：partial 不写 complete；only_qc 不引用未运行模块；
artifact_index 为显式文件集合（禁止递归 glob）。

错误码：0=成功；2=格式/配置/输入错误。
"""

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

EXIT_OK = 0
EXIT_INPUT_ERROR = 2


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="export_project_summary.py",
        description="项目级摘要与显式产物索引导出。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--effective-config", required=True, help="生效配置快照 YAML")
    parser.add_argument("--plan-json", required=True,
                        help="运行计划 JSON（module_status/analysis_targets/trait_metadata 路径）")
    parser.add_argument("--qc-audit", required=True, help="03.genotype_qc/qc_audit.json")
    parser.add_argument("--gate-json", required=True, help="01.gate/input_gate.json")
    parser.add_argument("--cohort-fam", required=True, help="最终 cohort.fam")
    parser.add_argument("--cohort-bim", required=True, help="最终 cohort.bim")
    parser.add_argument("--sample-order", required=True,
                        help="样本顺序 TSV（项目级权威：kinship sample_order）")
    parser.add_argument("--pipeline-version", required=True)
    parser.add_argument("--genome-version", required=True)
    parser.add_argument("--species", default="unspecified")
    parser.add_argument("--project-name", default="gwasflow_project")
    parser.add_argument("--seed", type=int, default=20261023)
    parser.add_argument("--only-qc", action="store_true",
                        help="only_qc 模式：不引用未运行的分析模块")
    parser.add_argument("--schema-version", default="gwasflow.project.v1")
    parser.add_argument("--out", required=True, help="project_summary.json 输出路径")
    parser.add_argument("--artifact-index-out", required=True, help="artifact_index.tsv 输出路径")
    return parser.parse_args(argv)


def fail(message, code=EXIT_INPUT_ERROR):
    print(f"[export_project_summary ERROR] {message}", file=sys.stderr)
    sys.exit(code)


def md5_of_file(path):
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_of_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path, label):
    p = Path(path)
    if not p.exists():
        fail(f"{label} 不存在: {p}")
    with open(p, encoding="utf-8") as handle:
        return json.load(handle)


def git_commit():
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True,
            cwd=str(Path(__file__).resolve().parent), check=False, timeout=10,
        )
        if result.returncode == 0:
            return result.stdout.strip()
        return "unavailable（git 仓库不可用）"
    except Exception:  # noqa: BLE001
        return "unavailable（git 调用失败）"


def count_lines(path):
    if not Path(path).exists():
        return None
    with open(path, encoding="utf-8") as handle:
        return sum(1 for line in handle if line.strip())


def main(argv=None):
    args = parse_args(argv)

    plan = load_json(args.plan_json, "运行计划 JSON")
    qc_audit = load_json(args.qc_audit, "qc_audit")
    gate = load_json(args.gate_json, "input_gate")
    module_status = plan.get("module_status", {})
    analysis_targets = plan.get("analysis_targets", [])
    trait_metadata_paths = plan.get("trait_metadata", {})

    # ---- 分析索引（每性状状态：从真实 metadata 读取，禁手工空 PASS）----
    analysis_index = []
    per_trait_counts = {}
    warnings = []
    for trait, meta_path in sorted(trait_metadata_paths.items()):
        meta = load_json(meta_path, f"性状 {trait} 的 harmonize metadata")
        status = meta.get("status", "unknown")
        counts = meta.get("counts", {})
        analysis_index.append({
            "trait": trait,
            "status": status,
            "n_valid_samples": None,   # 由准备 metadata 填充
            "m_tests": counts.get("n_valid_p"),
            "n_significant": counts.get("n_significant"),
            "n_unreported": counts.get("n_unreported"),
            "threshold": meta.get("multiple_testing", {}).get("threshold"),
            "metadata_path": meta_path,
        })
        per_trait_counts[trait] = {
            "m_tests": counts.get("n_valid_p"),
            "n_unreported": counts.get("n_unreported"),
        }
        if status == "partial":
            warnings.append(
                f"性状 {trait} 状态 partial（数值异常行见 invalid_variants.tsv）"
            )
        if status == "failed_software":
            warnings.append(
                f"性状 {trait} 上游软件失败（failed_software；退出码 "
                f"{meta.get('failure', {}).get('upstream_exit_code')}）"
            )

    # 准备 metadata（有效样本数）——trait_metadata 中可能附 prepare metadata 路径
    for entry in analysis_index:
        prep_path = plan.get("prepare_metadata", {}).get(entry["trait"])
        if prep_path and Path(prep_path).exists():
            prep_meta = load_json(prep_path, "prepare metadata")
            entry["n_valid_samples"] = prep_meta.get("n_valid")
            if prep_meta.get("pc_skipped"):
                warnings.append(
                    f"性状 {entry['trait']} PC 被跳过（only_qc 降级：pc_skipped）"
                )

    # ---- 模块运行状态（从产物存在性推断；only_qc 不引用未运行模块）----
    module_report = {}
    for module, info in module_status.items():
        enabled = info.get("enabled", False)
        module_report[module] = {
            "enabled": enabled,
            "auto_enable_reasons": info.get("auto_enable_reasons", []),
            "run_status": (
                "not_applicable" if not enabled else
                "executed" if analysis_targets else "executed"
            ),
        }

    # ---- 计数 ----
    n_samples_final = count_lines(args.cohort_fam)
    n_variants_final = count_lines(args.cohort_bim)
    stage_losses = [
        {
            "stage": s.get("stage"),
            "n_samples_removed": s.get("n_samples_removed"),
            "n_variants_removed": s.get("n_variants_removed"),
        }
        for s in qc_audit.get("stages", [])
    ]
    if qc_audit.get("warnings"):
        warnings.extend(qc_audit["warnings"])

    # ---- 七类审计信息 ----
    provenance = {
        "code_commit": git_commit(),
        "software_versions": {
            "gemma": "占位（M3：外部软件二进制未定，未验证项）",
            "plink2": "占位（M3：外部软件二进制未定，未验证项）",
            "bcftools": "占位（M3：外部软件二进制未定，未验证项）",
            "pipeline": args.pipeline_version,
        },
        "reference_versions": {
            "genome_version": args.genome_version,
            "note": "M2：生菜组装版本为占位声明，部署前须核对真实组装",
        },
        "seed": args.seed,
    }

    audit = {
        "input_digest": {
            "config_sha256": gate.get("config_sha256"),
            "inputs": [
                {"path": i.get("path"), "md5": i.get("md5")}
                for i in gate.get("inputs", [])
            ],
        },
        "config_snapshot": {
            "path": str(Path(args.effective_config).resolve()),
            "sha256": sha256_of_file(args.effective_config),
        },
        "code_commit": provenance["code_commit"],
        "software_environment": provenance["software_versions"],
        "reference_version": provenance["reference_versions"],
        "sample_order": {
            "path": str(Path(args.sample_order).resolve()),
            "sha256": sha256_of_file(args.sample_order) if Path(args.sample_order).exists() else "unavailable",
            "n_samples": count_lines(args.sample_order),
        },
    }

    # ---- artifact_index.tsv（显式文件集合，禁止递归 glob）----
    artifact_rows = []
    missing_artifacts = []
    for rel_path in analysis_targets:
        if not Path(rel_path).exists():
            missing_artifacts.append(rel_path)
            continue
        artifact_rows.append({
            "path": rel_path,
            "role": _role_of(rel_path),
            "size_bytes": os.path.getsize(rel_path),
            "md5": md5_of_file(rel_path),
            "source": "analysis_target（TargetPlan 显式集合）",
        })
    if missing_artifacts:
        warnings.append(
            f"目标清单中 {len(missing_artifacts)} 个产物缺失（部分失败或未运行模块）"
        )

    with open(args.artifact_index_out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["path", "role", "size_bytes", "md5", "source"],
            delimiter="\t",
        )
        writer.writeheader()
        for row in artifact_rows:
            writer.writerow(row)

    # ---- 未验证项清单（M3–M6 降级记录；不冒称已验证）----
    unverified_items = [
        "M3: PLINK2/GEMMA/BCFtools 真实二进制数值对账未执行（版本占位）",
        "M4: 真实植物群体端到端验收未执行（合成数据驱动 DAG/行为验证）",
        "M5: rich-loguru 日志插件生产接线未验收（三级降级链兜底）",
        "M6: BioReport GWAS 渲染器未验收（report=false；结构化摘要+静态图表）",
    ]
    if args.only_qc:
        unverified_items.append("only_qc 模式：关联分析未运行（DAG 截断）")

    # ---- 状态汇总（partial 不写 complete）----
    overall_status = "complete"
    if any(e["status"] in ("partial", "failed_software") for e in analysis_index):
        overall_status = "partial"
    if missing_artifacts:
        overall_status = "partial"

    summary = {
        "schema": args.schema_version,
        "project": {
            "project_name": args.project_name,
            "species": args.species,
            "genome_version": args.genome_version,
            "pipeline_version": args.pipeline_version,
            "only_qc": bool(args.only_qc),
        },
        "overall_status": overall_status,
        "module_status": module_report,
        "analysis_index": analysis_index,
        "counts": {
            "n_samples_final": n_samples_final,
            "n_variants_final": n_variants_final,
            "stage_losses": stage_losses,
            "per_trait_m_tests": {t: c["m_tests"] for t, c in per_trait_counts.items()},
            "per_trait_n_unreported": {t: c["n_unreported"] for t, c in per_trait_counts.items()},
        },
        "parameters_effective": {
            "effective_config": str(Path(args.effective_config).resolve()),
            "config_sha256": audit["config_snapshot"]["sha256"],
        },
        "provenance": provenance,
        "audit": audit,
        "warnings": warnings,
        "unverified_items": unverified_items,
        "artifact_index": {
            "path": str(Path(args.artifact_index_out).resolve()),
            "n_files": len(artifact_rows),
            "note": "显式文件集合（TargetPlan 派生），禁止递归 glob",
        },
    }
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)

    print(
        f"[export_project_summary] 摘要完成：整体状态 {overall_status}，"
        f"{len(artifact_rows)} 产物入索引，{len(warnings)} 警告"
    )
    return EXIT_OK


def _role_of(path):
    if path.startswith("01.gate/"):
        return "input_gate"
    if path.startswith("02.genotype/"):
        return "genotype_normalized"
    if path.startswith("03.genotype_qc/"):
        return "genotype_qc"
    if path.startswith("04.structure/"):
        return "population_structure"
    if path.startswith("05.phenotype/"):
        return "trait_preparation"
    if path.startswith("06.kinship/"):
        return "kinship"
    if path.startswith("07.association/"):
        return "association_raw"
    if path.startswith("08.results/"):
        return "unified_statistics"
    if path.startswith("09.loci/"):
        return "loci"
    if path.startswith("10.local_ld/"):
        return "regional_ld"
    if path.startswith("11.candidates/"):
        return "candidate_genes"
    if path.startswith("12.plots/"):
        return "visualization"
    if path.startswith("13.summary/"):
        return "project_summary"
    return "other"


if __name__ == "__main__":
    sys.exit(main())
