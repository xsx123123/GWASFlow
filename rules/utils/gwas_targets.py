#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GWASFlow TargetPlan 纯函数目标生成系统。

契约（spec §5.12.1 / 框架规范 §5.6）：
- 纯函数：不改写 config、不做 IO、不 sleep（可重复调用结果一致）
- 唯一调用点：snakefile 层 ALL_TARGETS = DataDeliver(...) 一次计算
- 15 个模块开关显式默认值 + only_qc 截断 + 依赖闭包不动点补齐
- plan_all 不含 deliver/report 规则自身输出以外的消费者输出（无循环依赖）
"""

from typing import NamedTuple

# 15 个模块开关（与 config.yaml / schema 逐字一致）
MODULE_SWITCHES = [
    "input_gate",
    "genotype_normalize",
    "genotype_qc",
    "population_structure",
    "phenotype_qc",
    "kinship_analysis",
    "association_analysis",
    "result_harmonization",
    "locus_analysis",
    "local_ld",
    "candidate_genes",
    "visualization",
    "export_summary",
    "deliver",
    "report",
]

# 基础模块（不可关；schema const: true）
BASE_MODULES = {"input_gate", "genotype_normalize", "genotype_qc", "export_summary"}

# only_qc 截断集合（强制关闭）
ONLY_QC_CUTOFF = [
    "population_structure",
    "kinship_analysis",
    "association_analysis",
    "result_harmonization",
    "locus_analysis",
    "local_ld",
    "candidate_genes",
    "visualization",
    "report",
]

# 模块依赖表（design §2.1.2-D；不含 PC 条件依赖，由 use_pcs 处理）
MODULE_DEPENDENCIES = {
    "input_gate": [],
    "genotype_normalize": ["input_gate"],
    "genotype_qc": ["genotype_normalize"],
    "population_structure": ["genotype_qc"],
    "phenotype_qc": ["genotype_qc"],
    "kinship_analysis": ["genotype_qc"],
    "association_analysis": ["genotype_qc", "phenotype_qc", "kinship_analysis"],
    "result_harmonization": ["association_analysis"],
    "locus_analysis": ["result_harmonization"],
    "local_ld": ["locus_analysis"],
    "candidate_genes": ["locus_analysis"],
    "visualization": ["result_harmonization", "phenotype_qc"],
    "export_summary": [],
    "deliver": ["export_summary"],
    "report": ["export_summary"],
}


class ModuleStatus(NamedTuple):
    module: str
    enabled: bool
    auto_enable_reasons: tuple  # tuple[str, ...]


class TargetPlan(NamedTuple):
    plan_analysis: tuple      # 12 个分析模块输出路径（启用者）
    plan_summary: tuple       # project_summary + artifact_index + effective_config
    plan_report: tuple        # 报告产物（report=true 时）
    plan_delivery: tuple      # delivery manifest 三件套（deliver=true 时）
    plan_all: tuple           # 上述四类并集（rule all 的 input）
    module_status: dict       # dict[str, ModuleStatus]


def resolve_modules(config, use_pcs=None):
    """解析生效模块开关（纯函数，只读 config）。

    Args:
      config: 深合并后的配置字典（只读，不修改）
      use_pcs: 是否有性状声明 n_pcs>0（True 时 phenotype_qc 需 population_structure；
               None 时跳过该条件依赖判断，用于 Step 6 开关级预览）

    Returns:
      dict[str, ModuleStatus]
    """
    status = {}
    for module in MODULE_SWITCHES:
        raw = config.get(module)
        if module in BASE_MODULES:
            enabled = True  # 基础模块不可关（schema 层已 const: true）
            if raw is False:
                # 防御：analysisyaml 违规覆盖时仍按基础模块处理并记录
                status[module] = ModuleStatus(
                    module, True, ("基础模块不可关闭（已忽略显式 false）",)
                )
                continue
        else:
            enabled = bool(raw)
        status[module] = ModuleStatus(module, enabled, ())

    # only_qc 截断
    if config.get("only_qc"):
        for module in ONLY_QC_CUTOFF:
            if status[module].enabled:
                status[module] = ModuleStatus(
                    module,
                    False,
                    status[module].auto_enable_reasons
                    + ("only_qc 模式截断（仅保留 gate/normalize/qc/phenotype_qc/summary/deliver）",),
                )

    # PC 条件依赖：有性状声明 n_pcs>0 且 phenotype_qc 启用时需要 population_structure
    # only_qc 模式下不自动启用（截断优先；n_pcs>0 的性状走 --skip-pcs 降级，D3）
    if (use_pcs and not config.get("only_qc")
            and status["phenotype_qc"].enabled
            and not status["population_structure"].enabled):
        status["population_structure"] = ModuleStatus(
            "population_structure",
            True,
            status["population_structure"].auto_enable_reasons
            + ("phenotype_qc 声明 n_pcs>0 需要 PCA 产物（自动启用）",),
        )

    # 依赖闭包（while 不动点补齐缺失上游）
    changed = True
    while changed:
        changed = False
        for module in MODULE_SWITCHES:
            if not status[module].enabled:
                continue
            for upstream in MODULE_DEPENDENCIES.get(module, []):
                if not status[upstream].enabled:
                    status[upstream] = ModuleStatus(
                        upstream,
                        True,
                        status[upstream].auto_enable_reasons
                        + (f"下游模块 {module} 启用需要 {upstream}（自动启用）",),
                    )
                    changed = True

    return status


def _collect_plan(config, samples, traits, module_status):
    """按模块收集器计算五类目标路径（纯函数）。"""
    from rules.utils.datadeliver import collect_all_targets

    return collect_all_targets(config, samples, traits, module_status)


def build_target_plan(config, samples, traits):
    """构建不可变 TargetPlan（纯函数，一次计算）。

    前置：config 通过 schema 与开关校验；samples/traits 已加载。
    """
    use_pcs = any(
        (t.get("n_pcs") or 0) > 0 for t in (traits or {}).values()
    )
    module_status = resolve_modules(config, use_pcs=use_pcs)

    analysis, summary, report, delivery = _collect_plan(
        config, samples, traits, module_status
    )

    all_targets = list(analysis) + list(summary) + list(report) + list(delivery)

    return TargetPlan(
        plan_analysis=tuple(analysis),
        plan_summary=tuple(summary),
        plan_report=tuple(report),
        plan_delivery=tuple(delivery),
        plan_all=tuple(all_targets),
        module_status=module_status,
    )


def format_module_status(module_status):
    """输出生效开关清单与自动启用原因（控制台/日志，spec §4.4-3/4）。"""
    lines = []
    for module in MODULE_SWITCHES:
        entry = module_status[module]
        state = "enabled" if entry.enabled else "disabled"
        line = f"  {module:<24} {state}"
        if entry.auto_enable_reasons:
            line += f"  <- {'; '.join(entry.auto_enable_reasons)}"
        lines.append(line)
    return "\n".join(lines)
