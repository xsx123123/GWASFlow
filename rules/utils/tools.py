#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GWASFlow 命令构建与执行辅助（CLI 横切约定）。

- build_command: 命令数组构建（禁止自由文本拼接，实现文档 §15）
- thread_env: BLAS/OpenMP 受控线程环境变量（spec §4.1-2）
- run_logged: 子进程执行，规则 log 与工具 --log-file 分离、退出码透传
- 错误码常量: 0=成功 2=格式/配置/输入错误 3=已导出结果但含数值异常
- ensure_empty_dir: 输出目录"新目录或空目录"保护
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

# CLI 横切约定：统一错误码
EXIT_OK = 0
EXIT_INPUT_ERROR = 2
EXIT_PARTIAL = 3


def build_command(template, **params):
    """命令数组构建：template 为字符串列表，占位符 {name} 由 params 填充。

    禁止自由文本拼接：所有参数必须经本函数显式注入。
    未知占位符或缺失参数立即报错（防参数静默丢失）。
    """
    command = []
    for item in template:
        if isinstance(item, str) and item.startswith("{") and item.endswith("}"):
            key = item[1:-1]
            if key not in params:
                raise KeyError(f"build_command 缺少参数 '{key}'（模板: {template}）")
            value = params[key]
            if value is None:
                continue
            command.append(str(value))
        else:
            command.append(str(item))
    return command


def thread_env(threads):
    """生成 BLAS/OpenMP 受控线程环境变量（与资源档位 threads 一致）。"""
    return {
        "OMP_NUM_THREADS": str(threads),
        "OPENBLAS_NUM_THREADS": str(threads),
        "MKL_NUM_THREADS": str(threads),
    }


def run_logged(argv, log, tool_log=None, extra_env=None, cwd=None):
    """执行子进程命令数组。

    - 规则 log 与工具自身 --log-file 必须为不同路径（双写交错防护）
    - stdout/stderr 追加写入规则 log
    - 退出码透传（外部软件非零退出不吞掉）
    """
    log_path = Path(log)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    if extra_env:
        env.update({k: str(v) for k, v in extra_env.items()})

    with open(log_path, "a", encoding="utf-8") as handle:
        handle.write(f"\n=== RUN: {' '.join(map(str, argv))} ===\n")
        handle.flush()
        result = subprocess.run(
            list(map(str, argv)),
            stdout=handle,
            stderr=subprocess.STDOUT,
            env=env,
            cwd=cwd,
            check=False,
        )
    return result.returncode


def ensure_empty_dir(out_dir):
    """输出目录必须为新目录或空目录（已有输出保护）。

    返回规范化绝对路径字符串；非空目录退出 2。
    """
    path = Path(out_dir).resolve()
    if path.exists():
        if not path.is_dir():
            print(
                f"[tools ERROR] 输出路径已存在且不是目录: {path}",
                file=sys.stderr,
            )
            sys.exit(EXIT_INPUT_ERROR)
        if any(path.iterdir()):
            print(
                f"[tools ERROR] 输出目录非空（拒绝覆盖既有产物，"
                f"重试请使用新的临时目录）: {path}",
                file=sys.stderr,
            )
            sys.exit(EXIT_INPUT_ERROR)
    else:
        path.mkdir(parents=True, exist_ok=True)
    return str(path)


def atomic_commit(tmp_dir, final_dir):
    """临时目录全部产物生成且校验通过后原子移动到声明输出目录。"""
    tmp = Path(tmp_dir)
    final = Path(final_dir)
    final.parent.mkdir(parents=True, exist_ok=True)
    if final.exists():
        if any(final.iterdir()):
            raise RuntimeError(f"目标目录非空，拒绝原子移动: {final}")
        final.rmdir()
    shutil.move(str(tmp), str(final))
    return str(final)


def sha256_of_file(path):
    """文件 SHA-256 摘要（审计与指纹绑定）。"""
    import hashlib

    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def md5_of_file(path):
    """文件 MD5 摘要（门控与交付校验）。"""
    import hashlib

    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
