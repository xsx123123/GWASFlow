#!/usr/bin/env python3
"""GWASFlow 资源档位管理（复用 FlowFrame 模板逻辑，去除插件依赖）。

接口契约: rule_resource(config, profile_name) -> {threads, mem_mb, runtime, queue?}
- 档位缺失 ValueError（带可用档位列表）
- local 模式 + skip_queue_on_local=True 时剥离 queue 参数
- 返回档位副本（无副作用）
- 线程申请不超过档位上限
"""


def rule_resource(config, profile_name, skip_queue_on_local=False, logger=None):
    """根据 config 获取规则的资源配置。"""
    queue_name = config.get("queue_id")

    cluster_config = config.get("cluster_config")
    if not cluster_config:
        error_msg = (
            "Critical Config Error: 'cluster_config' section is missing. "
            "Please check your configuration file."
        )
        if logger:
            logger.error(error_msg)
        raise ValueError(error_msg)

    resource_profiles = cluster_config.get("resource_profiles")
    if not resource_profiles:
        error_msg = (
            "Critical Config Error: 'resource_profiles' is missing within "
            "'cluster_config'. Please define resource profiles."
        )
        if logger:
            logger.error(error_msg)
        raise ValueError(error_msg)

    if profile_name not in resource_profiles:
        available = list(resource_profiles.keys())
        error_msg = (
            f"Resource profile '{profile_name}' not found in config. "
            f"Available: {available}"
        )
        if logger:
            logger.error(error_msg)
        raise ValueError(error_msg)

    profile = resource_profiles[profile_name].copy()

    execution_mode = config.get("execution_mode")
    current_cluster_name = cluster_config.get("queue_id")

    is_local_execution = (execution_mode == "local") or (
        current_cluster_name == "default"
    )

    if is_local_execution and skip_queue_on_local:
        if queue_name:
            if logger:
                logger.warning(
                    f"[Config Conflict] Execution mode is 'local', but a queue "
                    f"override ('{queue_name}') was provided. Ignored queue "
                    f"parameter to ensure local execution success."
                )
        else:
            if logger:
                logger.debug(
                    f"Running locally. Stripping queue params for '{profile_name}'."
                )
        profile.pop("queue", None)
        profile.pop("queue_type", None)
        return profile

    if queue_name:
        profile["queue"] = queue_name
    elif "queue_type" in profile:
        clusters = cluster_config.get("clusters", {})
        current_cluster = clusters.get(
            current_cluster_name, clusters.get("default", {})
        )
        queues = current_cluster.get("queues", {})
        queue_key = profile["queue_type"]
        if queue_key in queues:
            profile["queue"] = queues[queue_key]
            if logger:
                logger.debug(
                    f"Mapped queue_type '{queue_key}' -> '{profile['queue']}'"
                )
        else:
            if logger:
                logger.warning(
                    f"Queue type '{queue_key}' not defined in cluster "
                    f"'{current_cluster_name}'. Job submits without explicit queue."
                )

    return profile
