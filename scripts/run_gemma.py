#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GWASFlow GEMMA 关联分析执行器（参数手册 §4 契约生产化，D4/D5/D10）。

命令数组构建（内部参数映射表全部落实）：
  -bfile cohort 全集（D10：与独立对照命令同参数可比；MAC 差异经 unreported 对账暴露）
  -p 完整 N 行外部表型（含 NA 掩码，不读 FAM 原有表型）
  -c 显式截距协变量（covariates.gemma.txt 已含截距）
  -k 全队列 K
  -lmm wald|lrt|score|all → 1|2|3|4（禁止逐位点取多检验最小 P）
  -maf/-miss 显式生效值；-hwe 0 显式关闭；-r2 0.9999 记录
  -outdir/-o 避免共享默认 output

红线：
- K 顺序与性状 sample_order 对账失败退出 2
- BLAS/OpenMP 线程环境与资源档位 threads 一致（thread_env 注入）
- GEMMA 非零退出透传退出码，日志保留原始错误
- run_metadata 绑定五类信息（BED/BIM/FAM、trait metadata、K、assoc、二进制摘要）

错误码：0=成功；2=格式/配置/输入错误；GEMMA 退出码透传。
"""

import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

EXIT_OK = 0
EXIT_INPUT_ERROR = 2
LMM_MAP = {"wald": "1", "lrt": "2", "score": "3", "all": "4"}


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="run_gemma.py",
        description="GEMMA LMM 关联分析（每性状独立运行；性状间无依赖）。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--bfile", required=True, help="最终 cohort BED 前缀（全集）")
    parser.add_argument("--trait-dir", required=True, help="05.phenotype/{trait}/ 准备产物目录")
    parser.add_argument("--trait", required=True, help="性状 ID")
    parser.add_argument("--kinship", required=True, help="06.kinship/kinship.cXX.txt（或 sXX）")
    parser.add_argument("--kinship-sample-order", required=True,
                        help="06.kinship/sample_order.tsv（K 顺序对账）")
    parser.add_argument("--kinship-type", default="centered",
                        choices=["centered", "standardized"], help="K 估计方法")
    parser.add_argument("--trait-mac", required=True, help="07.association/{trait}/trait_mac.tsv")
    parser.add_argument("--test", default="wald", choices=["wald", "lrt", "score", "all"],
                        help="主检验（预声明；all 时各检验独立归档，禁止逐位点取最小 P）")
    parser.add_argument("--maf", type=float, default=0.05, help="GEMMA -maf 显式生效值（比例）")
    parser.add_argument("--miss", type=float, default=0.10, help="GEMMA -miss 阈值（比例）")
    parser.add_argument("--threads", type=int, default=1,
                        help="BLAS/OpenMP 线程上限（注入环境变量；非 GEMMA 原生线程）")
    parser.add_argument("--resource-profile", default="medium_resource",
                        help="实际生效资源档位名（记入 run_metadata）")
    parser.add_argument("--keep-going", action="store_true",
                        help="规则层失败传导：GEMMA 非零退出时补全空表头后退出 0"
                             "（真实退出码记录于 run_metadata；CLI 独立调用默认透传）")
    parser.add_argument("--out-dir", required=True, help="输出目录（必须为新目录或空目录）")
    return parser.parse_args(argv)


def fail(message, code=EXIT_INPUT_ERROR):
    print(f"[run_gemma ERROR] {message}", file=sys.stderr)
    sys.exit(code)


def ensure_empty_dir(out_dir):
    path = Path(out_dir).resolve()
    if path.exists():
        if not path.is_dir() or any(path.iterdir()):
            fail(f"输出目录非空或非目录（拒绝覆盖既有产物）: {path}")
    else:
        path.mkdir(parents=True, exist_ok=True)
    return path


def sha256_of_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_order_tsv(path):
    entries = []
    with open(path, encoding="utf-8") as handle:
        header = handle.readline()
        if not header.lower().startswith("fid"):
            fail(f"顺序文件表头异常: {path}")
        for line in handle:
            parts = line.split()
            if len(parts) >= 2:
                entries.append((parts[0], parts[1]))
    return entries


def count_file_lines(path):
    with open(path, encoding="utf-8") as handle:
        return sum(1 for line in handle if line.strip())


def main(argv=None):
    args = parse_args(argv)
    out_dir = ensure_empty_dir(args.out_dir)
    trait_dir = Path(args.trait_dir)

    pheno = trait_dir / "phenotype.gemma.txt"
    covar = trait_dir / "covariates.gemma.txt"
    trait_meta_path = trait_dir / "metadata.json"
    trait_order = trait_dir / "sample_order.tsv"

    for path in (pheno, covar, trait_meta_path, trait_order, Path(args.kinship),
                 Path(args.kinship_sample_order), Path(args.trait_mac),
                 Path(f"{args.bfile}.bed"), Path(f"{args.bfile}.bim"),
                 Path(f"{args.bfile}.fam")):
        if not path.exists():
            fail(f"必需输入不存在: {path}")

    # ---- K 顺序与性状 sample_order 对账 ----
    kinship_order = read_order_tsv(args.kinship_sample_order)
    trait_order_entries = read_order_tsv(trait_order)
    if kinship_order != trait_order_entries:
        fail(
            f"K 的样本顺序与性状 sample_order 不一致（对账失败）。"
            f"K 首行 {kinship_order[:2]} vs 性状首行 {trait_order_entries[:2]}；"
            f"长度 {len(kinship_order)} vs {len(trait_order_entries)}"
        )

    # ---- 表型/协变量行数一致性（全矩阵路线）----
    n_fam = count_file_lines(f"{args.bfile}.fam")
    n_pheno = count_file_lines(pheno)
    n_covar = count_file_lines(covar)
    if not (n_pheno == n_covar == n_fam):
        fail(
            f"表型/协变量行数与 FAM 不一致: pheno={n_pheno}, covar={n_covar}, "
            f"fam={n_fam}（全矩阵+缺失掩码要求完整 N 行）"
        )

    # ---- 命令数组构建（线程环境受控）----
    env = os.environ.copy()
    thread_envs = {
        "OMP_NUM_THREADS": str(args.threads),
        "OPENBLAS_NUM_THREADS": str(args.threads),
        "MKL_NUM_THREADS": str(args.threads),
    }
    env.update(thread_envs)

    gemma_prefix = out_dir / "association"
    argv_cmd = [
        "gemma",
        "-bfile", args.bfile,
        "-p", str(pheno),
        "-c", str(covar),
        "-k", args.kinship,
        "-lmm", LMM_MAP[args.test],
        "-maf", str(args.maf),
        "-miss", str(args.miss),
        "-hwe", "0",
        "-r2", "0.9999",
        "-n", "1",
        "-outdir", str(out_dir),
        "-o", gemma_prefix.name,
    ]

    log_path = out_dir / "gemma.run.log"
    with open(log_path, "w", encoding="utf-8") as handle:
        handle.write(f"=== RUN: {' '.join(argv_cmd)} ===\n")
        handle.write(f"=== thread_env: {thread_envs} ===\n")
        handle.flush()
        result = subprocess.run(argv_cmd, stdout=handle, stderr=subprocess.STDOUT,
                                env=env, check=False)
    exit_code = result.returncode

    # ---- GEMMA 自身 log 副本（失败时保留空文件以维持 output 契约）----
    gemma_log_src = out_dir / f"{gemma_prefix.name}.log.txt"
    if gemma_log_src.exists():
        os.replace(gemma_log_src, out_dir / "association.log.txt")
    else:
        (out_dir / "association.log.txt").touch()

    # ---- gemma.help.txt（二进制摘要：无参数运行输出版本与用法）----
    help_log = out_dir / "gemma_help_capture.log"
    with open(help_log, "w", encoding="utf-8") as handle:
        help_result = subprocess.run(["gemma"], stdout=handle,
                                     stderr=subprocess.STDOUT, check=False)
    with open(help_log, encoding="utf-8", errors="replace") as handle:
        help_text = handle.read()
    with open(out_dir / "gemma.help.txt", "w", encoding="utf-8") as handle:
        handle.write(help_text[:20000])
    os.remove(help_log)

    # ---- 原始输出归档（按检验类型重命名；all 时三件独立归档）----
    assoc_outputs = {}
    _empty_header = ("chr\trs\tn_miss\tn_obs\tallele1\tallele2\taf"
                     "\tbeta\tse\tlogl_H1\tl_remle\tp_wald\n")
    expected_files = {
        "wald": [f"{gemma_prefix.name}.assoc.txt"],
        "lrt": [f"{gemma_prefix.name}.lrt.txt"],
        "score": [f"{gemma_prefix.name}.score.txt"],
        "all": [
            f"{gemma_prefix.name}.assoc.txt",
            f"{gemma_prefix.name}.lrt.txt",
            f"{gemma_prefix.name}.score.txt",
        ],
    }[args.test]
    for src_name in expected_files:
        src = out_dir / src_name
        if src.exists():
            dst = out_dir / ("association" + src_name[len(gemma_prefix.name):])
            os.replace(src, dst)
            assoc_outputs[dst.name] = str(dst)
        elif exit_code != 0:
            # 失败传导：补全空表头（DAG 可恢复；其余性状不受影响）
            dst = out_dir / ("association" + src_name[len(gemma_prefix.name):])
            with open(dst, "w", encoding="utf-8") as handle:
                handle.write(_empty_header)
            assoc_outputs[dst.name] = str(dst)

    # ---- run_metadata.json（五类绑定 + 软件四过滤生效值 + 资源档位）----
    if exit_code != 0:
        status = "failed_software"
    else:
        status = "software_exited_zero"

    with open(trait_meta_path, encoding="utf-8") as handle:
        trait_metadata = json.load(handle)

    run_metadata = {
        "schema": "gwasflow.gemma_run.v1",
        "trait": args.trait,
        "test": args.test,
        "lmm_code": LMM_MAP[args.test],
        "status": status,
        "exit_code": exit_code,
        "bindings": {
            "bed": {"path": str(Path(f"{args.bfile}.bed").resolve()),
                    "sha256": sha256_of_file(f"{args.bfile}.bed")},
            "bim": {"path": str(Path(f"{args.bfile}.bim").resolve()),
                    "sha256": sha256_of_file(f"{args.bfile}.bim")},
            "fam": {"path": str(Path(f"{args.bfile}.fam").resolve()),
                    "sha256": sha256_of_file(f"{args.bfile}.fam")},
            "trait_metadata": {
                "path": str(trait_meta_path.resolve()),
                "sha256": sha256_of_file(trait_meta_path),
                "n_valid": trait_metadata.get("n_valid"),
                "n_total": trait_metadata.get("n_total"),
                "trait": trait_metadata.get("trait"),
            },
            "kinship": {
                "path": str(Path(args.kinship).resolve()),
                "sha256": sha256_of_file(args.kinship),
                "type": args.kinship_type,
                "n_samples": len(kinship_order),
            },
            "assoc": {name: {"path": p, "sha256": sha256_of_file(p)}
                      for name, p in assoc_outputs.items()} or "no_output",
            "binary": {
                "name": "gemma",
                "help_digest": hashlib.sha256(help_text.encode()).hexdigest(),
                "help_file": "gemma.help.txt",
            },
            "trait_mac": {"path": str(Path(args.trait_mac).resolve()),
                          "sha256": sha256_of_file(args.trait_mac)},
        },
        "software_filters": {
            "maf": args.maf,
            "miss": args.miss,
            "hwe": 0,
            "r2": 0.9999,
            "note": "GEMMA 四过滤显式生效值（-hwe 0=关闭；-r2 为协变量相关过滤阈值，"
                    "非 LD clumping 阈值）",
        },
        "resource": {
            "profile": args.resource_profile,
            "threads": args.threads,
            "thread_env": thread_envs,
        },
        "command": " ".join(argv_cmd),
        "pheno_rows": n_pheno,
        "covar_rows": n_covar,
        "fam_rows": n_fam,
        "masking_route": "全队列 K + 性状缺失掩码（-p 完整 N 行含 NA）",
    }
    with open(out_dir / "run_metadata.json", "w", encoding="utf-8") as handle:
        json.dump(run_metadata, handle, indent=2, ensure_ascii=False)

    if exit_code != 0:
        message = (
            f"[run_gemma ERROR] GEMMA 退出码 {exit_code}；"
            f"日志保留原始错误: {log_path}"
        )
        print(message, file=sys.stderr)
        if args.keep_going:
            print(
                "[run_gemma RULE] --keep-going 失败传导：已补全空表头产物，"
                "状态 failed_software 透传至 run_metadata（供 project_summary 记录）",
                file=sys.stderr,
            )
            return EXIT_OK
        sys.exit(exit_code)

    print(
        f"[run_gemma] 性状 '{args.trait}' 关联分析完成（test={args.test}），"
        f"产物目录 {out_dir}"
    )
    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
