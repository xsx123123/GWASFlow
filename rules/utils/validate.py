#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GWASFlow 配置校验：analysisyaml 深合并 / Genome_Version / 未实现模块开关。

未实现功能（loco/farmcpu/blink/imputation/multi_environment/human_profile/
pgen_import/hapmap_import）启用即报错，错误信息指明"属后续工作包"（错误码 2）。
"""

import os
import sys

NOT_IMPLEMENTED_KEYS = [
    "loco",
    "farmcpu",
    "blink",
    "imputation",
    "multi_environment",
    "human_profile",
    "pgen_import",
    "hapmap_import",
]


def _die(message: str, code: int = 2) -> None:
    print(f"[validate ERROR] {message}", file=sys.stderr)
    sys.exit(code)


def load_user_config(config, cmd_arg_name="analysisyaml"):
    """解析 CLI 传入的 analysisyaml 并深合并进 config（最高优先级）。

    复用 snakemake.utils.update_config 实现递归合并；文件不存在退出 2。
    """
    from snakemake.utils import update_config

    custom_path = config.get(cmd_arg_name)
    if not custom_path:
        return

    if not os.path.exists(custom_path):
        _die(f"找不到指定的项目配置文件: {custom_path}，请检查路径。")

    import yaml

    print(f"[Config Info] 加载外部项目配置: {custom_path}")
    try:
        with open(custom_path, "r", encoding="utf-8") as handle:
            custom_data = yaml.safe_load(handle)
    except Exception as exc:  # noqa: BLE001
        _die(f"解析 YAML 失败: {exc}")

    if custom_data:
        update_config(config, custom_data)
    else:
        print(f"[Config Warning] 文件 {custom_path} 为空，跳过加载。")


def validate_genome_version(config, logger=None):
    """校验 Genome_Version 属于 can_use_genome_version 列表。"""
    if not config:
        raise ValueError("Config 字典不能为空。")

    raw_version = config.get("Genome_Version")
    allowed_list = config.get("can_use_genome_version")

    if not raw_version or not allowed_list:
        msg = "Config 缺少必需键: 'Genome_Version' 或 'can_use_genome_version'。"
        if logger:
            logger.error(msg)
        raise ValueError(msg)

    if not isinstance(allowed_list, list):
        msg = f"'can_use_genome_version' 必须为列表，实际为 {type(allowed_list)}。"
        if logger:
            logger.error(msg)
        raise TypeError(msg)

    clean_version = str(raw_version).strip()
    if clean_version in allowed_list:
        if logger:
            logger.info(f"Genome version verified: '{clean_version}' is supported.")
        return clean_version

    msg = (
        f"不支持的基因组版本: '{clean_version}'。"
        f"支持的版本: {allowed_list}"
    )
    if logger:
        logger.error(msg)
    raise ValueError(msg)


def validate_module_switches(config, logger=None):
    """未实现功能键值非 false/null 时报错（D8 第二层运行前拦截）。

    错误信息指明键名与"属后续工作包"，错误码 2。
    """
    violations = []
    for key in NOT_IMPLEMENTED_KEYS:
        value = config.get(key)
        if value is not None and value is not False:
            violations.append((key, value))

    if violations:
        lines = [
            f"  - {key} = {value!r}（该能力属后续工作包，本轮未实现，"
            f"请保持 false 或删除该键）"
            for key, value in violations
        ]
        msg = "检测到未实现功能被启用:\n" + "\n".join(lines)
        if logger:
            logger.error(msg)
        _die(msg)
