#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# GWASFlow - 01.common.smk: Common imports + 3-level logger fallback chain (no business rules)
#
# Contracts:
# - Logger 3-level fallback chain: rich-loguru plugin -> loguru -> stdlib logging (D9)
# - No top-level bare plugin imports (try/except fallback is a development safety net)
# - The logs/<NN>.<module>/ numbered-directory convention is established in this file
# - Production environments must install snakemake-logger-plugin-rich-loguru (M5 fallback note in README)

import os

# ---- Logger 3-level fallback chain ----
try:
    from snakemake_logger_plugin_rich_loguru import get_analysis_logger
    logger = get_analysis_logger()
    LOGGER_TYPE = "rich-loguru plugin"
except ImportError:
    try:
        from loguru import logger
        LOGGER_TYPE = "standard loguru"
        logger.warning(
            "rich-loguru plugin not installed, falling back to standard loguru. "
            "Production environments must install snakemake-logger-plugin-rich-loguru."
        )
    except ImportError:
        import logging

        logging.basicConfig(
            level=os.environ.get("GWASFLOW_LOG_LEVEL", "INFO")
        )
        logger = logging.getLogger("GWASFlow")
        LOGGER_TYPE = "stdlib logging"
        logger.warning(
            "Neither the rich-loguru plugin nor loguru is installed, falling back to stdlib logging. "
            "Production environments must install snakemake-logger-plugin-rich-loguru."
        )

logger.info(f"GWASFlow logger initialized ({LOGGER_TYPE})")

# ---- Common utility imports (reused by subsequent rule files) ----
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
