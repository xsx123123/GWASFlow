#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GWASFlow 性状级精确 MAC 计算（spec §5.7.1-4，D4）。

plink2 --keep keep.tsv --freq 在该性状有效材料上计算逐标记精确 allele counts，
生成允许测试标记集（MAC >= mac-min）。队列级 MAC 不能当成性状级精确 MAC。

错误码：0=成功；2=格式/配置/输入错误；plink2 非零退出透传。
"""

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

EXIT_OK = 0
EXIT_INPUT_ERROR = 2


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="compute_trait_mac.py",
        description="在该性状有效材料（keep.tsv）上计算逐标记精确 MAC 与允许测试标记集。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--bfile", required=True, help="最终 cohort BED 前缀")
    parser.add_argument("--fam", required=True, help="cohort.fam（身份对账基准）")
    parser.add_argument("--keep", required=True, help="该性状 keep.tsv（有效 FID/IID）")
    parser.add_argument("--mac-min", type=int, default=5, help="最小 MAC 阈值（正整数）")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--out-dir", required=True, help="输出目录（新目录或空目录）")
    return parser.parse_args(argv)


def fail(message, code=EXIT_INPUT_ERROR):
    print(f"[compute_trait_mac ERROR] {message}", file=sys.stderr)
    sys.exit(code)


def ensure_empty_dir(out_dir):
    path = Path(out_dir).resolve()
    if path.exists():
        if not path.is_dir() or any(path.iterdir()):
            fail(f"输出目录非空或非目录（拒绝覆盖既有产物）: {path}")
    else:
        path.mkdir(parents=True, exist_ok=True)
    return path


def read_fam_ids(fam_path):
    pairs = []
    with open(fam_path, encoding="utf-8") as handle:
        for line in handle:
            parts = line.split()
            if len(parts) >= 2:
                pairs.append((parts[0], parts[1]))
    return pairs


def main(argv=None):
    args = parse_args(argv)
    out_dir = ensure_empty_dir(args.out_dir)

    # keep 身份与 FAM 对账（不一致退出 2）
    fam_ids = set(read_fam_ids(args.fam))
    keep_ids = []
    with open(args.keep, encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            parts = line.split()
            if len(parts) < 2:
                fail(f"keep.tsv 第 {line_no} 行列数不足: {line.strip()!r}")
            keep_ids.append((parts[0], parts[1]))
    keep_set = set(keep_ids)
    if len(keep_set) != len(keep_ids):
        fail("keep.tsv 存在重复身份")
    unknown = keep_set - fam_ids
    if unknown:
        fail(f"keep.tsv 身份不在 FAM 中: {sorted(unknown)}")

    n_keep = len(keep_ids)
    log_path = out_dir / "plink2.log"

    # plink2 --keep --freq（精确 allele counts）
    freq_prefix = out_dir / "trait_freq"
    argv_cmd = [
        "plink2", "--bfile", args.bfile,
        "--keep", args.keep,
        "--freq",
        "--threads", str(args.threads),
        "--out", str(freq_prefix),
    ]
    with open(log_path, "w", encoding="utf-8") as handle:
        handle.write(f"=== RUN: {' '.join(argv_cmd)} ===\n")
        handle.flush()
        result = subprocess.run(argv_cmd, stdout=handle, stderr=subprocess.STDOUT,
                                check=False)
    if result.returncode != 0:
        fail(f"plink2 --keep --freq 失败（退出码 {result.returncode}），日志: {log_path}")

    afreq_path = freq_prefix.with_suffix(".afreq")
    if not afreq_path.exists():
        fail(f"plink2 未产出 .afreq 文件（期望 {afreq_path}）")

    # 解析 .afreq：#CHROM ID REF ALT ALT_FREQS [OBS_CT]
    rows = []
    with open(afreq_path, encoding="utf-8") as handle:
        header = [c.strip().lstrip("#") for c in handle.readline().split()]
        for line in handle:
            parts = line.strip().split()
            if not parts:
                continue
            record = dict(zip(header, parts))
            variant_id = record.get("ID", "")
            chrom = record.get("CHROM", "")
            ref, alt = record.get("REF", ""), record.get("ALT", "")
            try:
                alt_freq = float(record.get("ALT_FREQS", "nan"))
            except ValueError:
                alt_freq = float("nan")
            obs_ct = None
            if "OBS_CT" in record and record["OBS_CT"] not in (".", "NA"):
                try:
                    obs_ct = int(record["OBS_CT"])
                except ValueError:
                    obs_ct = None
            if obs_ct is None:
                obs_ct = 2 * n_keep
            alt_count = round(alt_freq * obs_ct)
            mac = min(alt_count, obs_ct - alt_count)
            pos = record.get("POS", "")
            rows.append({
                "variant_id": variant_id,
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "alt_freq": alt_freq,
                "alt_count": alt_count,
                "mac": mac,
                "allowed": mac >= args.mac_min,
            })

    with open(out_dir / "trait_mac.tsv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["variant_id", "chrom", "pos", "ref", "alt",
                         "alt_count", "mac", "allowed"])
        for row in rows:
            writer.writerow([
                row["variant_id"], row["chrom"], row["pos"], row["ref"], row["alt"],
                row["alt_count"], row["mac"],
                "true" if row["allowed"] else "false",
            ])

    n_allowed = sum(1 for r in rows if r["allowed"])
    summary = {
        "schema": "gwasflow.trait_mac.v1",
        "n_keep_samples": n_keep,
        "n_variants_total": len(rows),
        "n_allowed": n_allowed,
        "mac_min": args.mac_min,
        "method": "plink2 --keep <trait keep.tsv> --freq（精确 allele counts，D4）",
        "note": "队列级 MAC 不能当成性状级精确 MAC；GEMMA 将按其分析材料重新过滤",
    }
    with open(out_dir / "trait_mac_summary.json", "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)

    print(
        f"[compute_trait_mac] 精确 MAC 完成：{n_allowed}/{len(rows)} 标记允许测试"
        f"（mac_min={args.mac_min}，N_keep={n_keep}）"
    )
    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
