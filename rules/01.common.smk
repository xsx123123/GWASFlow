#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# GWASFlow - 01.common.smk: 公共导入 + logger 三级降级链（无业务 rule）
#
# 契约：
# - 日志三级降级链: rich-loguru 插件 -> loguru -> stdlib logging（D9）
# - 禁止顶层裸 import 插件（try/except fallback 为开发兜底）
# - logs/<NN>.<module>/ 编号目录约定在本文件定型
# - 生产环境必须安装 snakemake-logger-plugin-rich-loguru（M5 降级说明见 README）

import os

# ---- Logger 三级降级链 ----
try:
    from snakemake_logger_plugin_rich_loguru import get_analysis_logger
    logger = get_analysis_logger()
    LOGGER_TYPE = "rich-loguru plugin"
except ImportError:
    try:
        from loguru import logger
        LOGGER_TYPE = "standard loguru"
        logger.warning(
            "rich-loguru 插件未安装，降级到标准 loguru。"
            "生产环境必须安装 snakemake-logger-plugin-rich-loguru。"
        )
    except ImportError:
        import logging

        logging.basicConfig(
            level=os.environ.get("GWASFLOW_LOG_LEVEL", "INFO")
        )
        logger = logging.getLogger("GWASFlow")
        LOGGER_TYPE = "stdlib logging"
        logger.warning(
            "rich-loguru 插件与 loguru 均未安装，降级到标准库 logging。"
            "生产环境必须安装 snakemake-logger-plugin-rich-loguru。"
        )

logger.info(f"GWASFlow logger initialized ({LOGGER_TYPE})")

# ---- 公共工具导入（供后续规则文件复用） ----
from rules.utils.tools import (
    build_command,
    thread_env,
    run_logged,
    ensure_empty_dir,
    atomic_commit,
    sha256_of_file,
    md5_of_file,
    EXIT_OK,
    EXIT_INPUT_ERROR,
    EXIT_PARTIAL,
)
from rules.utils.resource_manager import rule_resource
