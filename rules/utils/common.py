#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GWASFlow 兼容入口：AnalysisTargets / DataDeliver 复用同一 TargetPlan。

框架规范 §5.6-2 唯一调用点约束：snakefile 层一次计算 ALL_TARGETS，
deliver/report 规则的 input 引用该计划字段，不再重算（规避模板已验证 bug）。
"""

from rules.utils.gwas_targets import build_target_plan


def AnalysisTargets(config, samples, traits):
    """返回分析类目标列表（plan_analysis）。"""
    plan = build_target_plan(config, samples, traits)
    return list(plan.plan_analysis)


def DataDeliver(config, samples, traits):
    """返回完整 TargetPlan（不可变五类目标）。"""
    return build_target_plan(config, samples, traits)


def ReportData(config, samples, traits):
    """返回报告类目标列表（report 开关由 plan 内部尊重）。"""
    plan = build_target_plan(config, samples, traits)
    return list(plan.plan_report)
