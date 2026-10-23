#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GWASFlow 全队列亲缘关系矩阵计算（参数手册 §9.3）。

GEMMA -gk 1|2 -p kinship_all.pheno.txt -bfile cohort：
- 全 1 表型强制覆盖全队列（不读 FAM 原有缺失表型）
- type=centered|standardized 分别对应 -gk 1|2（禁止跨项目/跨方法缓存混用）
- mode=genome_wide（loco 已在 Schema 层拒绝）
- 输出统一命名 + sample_order.tsv（与 FAM 顺序绑定）

错误码：0=成功；2=格式/配置/输入错误；GEMMA 非零退出透传。
"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

EXIT_OK = 0
EXIT_INPUT_ERROR = 2


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="compute_kinship.py",
        description="GEMMA 全队列 kinship 矩阵计算（genome_wide 模式）。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--bfile", required=True, help="最终 cohort BED 前缀")
    parser.add_argument("--pheno-all-ones", required=True,
                        help="全 1 表型文件（完整 N 行，强制全队列）")
    parser.add_argument("--fam", required=True, help="cohort.fam（样本顺序绑定）")
    parser.add_argument("--mode", default="genome_wide", choices=["genome_wide"],
                        help="模式（loco 属后续工作包，Schema 层拒绝）")
    parser.add_argument("--type", default="centered", choices=["centered", "standardized"],
                        help="K 估计方法（centered=-gk 1；standardized=-gk 2）")
    parser.add_argument("--threads", type=int, default=1,
                        help="BLAS/OpenMP 线程上限（注入环境变量）")
    parser.add_argument("--out-dir", required=True, help="输出目录（必须为新目录或空目录）")
    return parser.parse_args(argv)


def fail(message, code=EXIT_INPUT_ERROR):
    print(f"[compute_kinship ERROR] {message}", file=sys.stderr)
    sys.exit(code)


def ensure_empty_dir(out_dir):
    path = Path(out_dir).resolve()
    if path.exists():
        if not path.is_dir() or any(path.iterdir()):
            fail(f"输出目录非空或非目录（拒绝覆盖既有产物）: {path}")
    else:
        path.mkdir(parents=True, exist_ok=True)
    return path


def read_fam_order(fam_path):
    order = []
    with open(fam_path, encoding="utf-8") as handle:
        for line in handle:
            parts = line.split()
            if len(parts) >= 2:
                order.append((parts[0], parts[1]))
    return order


def main(argv=None):
    args = parse_args(argv)
    out_dir = ensure_empty_dir(args.out_dir)

    fam_order = read_fam_order(args.fam)
    n_expected = len(fam_order)
    with open(args.pheno_all_ones, encoding="utf-8") as handle:
        n_pheno = sum(1 for line in handle if line.strip())
    if n_pheno != n_expected:
        fail(
            f"全 1 表型文件行数 {n_pheno} != FAM N {n_expected}"
            f"（全队列 K 要求完整 N 行覆盖）"
        )

    gk = "1" if args.type == "centered" else "2"
    matrix_suffix = "cXX" if args.type == "centered" else "sXX"

    env = os.environ.copy()
    env.update({
        "OMP_NUM_THREADS": str(args.threads),
        "OPENBLAS_NUM_THREADS": str(args.threads),
        "MKL_NUM_THREADS": str(args.threads),
    })

    gemma_out_prefix = out_dir / "gemma_kinship"
    argv_cmd = [
        "gemma",
        "-bfile", args.bfile,
        "-gk", gk,
        "-p", args.pheno_all_ones,
        "-n", "1",
        "-outdir", str(out_dir),
        "-o", gemma_out_prefix.name,
    ]
    log_path = out_dir / "gemma.run.log"
    with open(log_path, "w", encoding="utf-8") as handle:
        handle.write(f"=== RUN: {' '.join(argv_cmd)} ===\n")
        handle.flush()
        result = subprocess.run(argv_cmd, stdout=handle, stderr=subprocess.STDOUT,
                                env=env, check=False)
    if result.returncode != 0:
        fail(
            f"GEMMA -gk {gk} 退出码 {result.returncode}（透传），"
            f"日志: {log_path}"
        )

    # GEMMA 输出重命名为统一名
    src_matrix = out_dir / f"{gemma_out_prefix.name}.{matrix_suffix}.txt"
    dst_matrix = out_dir / f"kinship.{matrix_suffix}.txt"
    if not src_matrix.exists():
        # 备选：GEMMA 有时直接输出 <prefix>.cXX.txt 于 outdir
        candidates = list(out_dir.glob(f"*.{matrix_suffix}.txt"))
        if not candidates:
            fail(f"GEMMA 未产出 K 矩阵（期望 {src_matrix}），日志: {log_path}")
        src_matrix = candidates[0]
    os.replace(src_matrix, dst_matrix)

    src_log = out_dir / f"{gemma_out_prefix.name}.log.txt"
    dst_log = out_dir / "kinship.log.txt"
    if src_log.exists():
        os.replace(src_log, dst_log)

    # sample_order.tsv（与 FAM 顺序绑定）
    with open(out_dir / "sample_order.tsv", "w", encoding="utf-8") as handle:
        handle.write("fid\tiid\n")
        for fid, iid in fam_order:
            handle.write(f"{fid}\t{iid}\n")

    # 运行侧信息
    side_info = {
        "schema": "gwasflow.kinship_run.v1",
        "mode": args.mode,
        "type": args.type,
        "gemma_gk": gk,
        "n_samples": n_expected,
        "bfile": str(Path(args.bfile).resolve()),
        "pheno_all_ones": str(Path(args.pheno_all_ones).resolve()),
        "threads_env": {"OMP_NUM_THREADS": str(args.threads),
                        "OPENBLAS_NUM_THREADS": str(args.threads),
                        "MKL_NUM_THREADS": str(args.threads)},
        "matrix": str(dst_matrix),
        "note": "全队列 K（全 1 表型强制覆盖）；缓存身份须含估计方法与样本集合",
    }
    with open(out_dir / "kinship_run_info.json", "w", encoding="utf-8") as handle:
        json.dump(side_info, handle, indent=2, ensure_ascii=False)

    print(f"[compute_kinship] 全队列 K 完成（type={args.type}，N={n_expected}）: {dst_matrix}")
    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
