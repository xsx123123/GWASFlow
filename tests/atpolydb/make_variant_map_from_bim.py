#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""AtPolyDB（TAIR9）测试数据辅助：从 TPED 或 BIM 生成 GWASFlow BED 路线 variant_map。

语义诚实声明（与数据包 README §4 一致）：
- PLINK BIM 的 A1/A2 是**观测等位基因**，无 REF/ALT 参考方向断言；
- 本工具按登记约定 ref=A2(allele2)、alt=A1(allele1)（plink --keep-allele-order
  保持 tped 列序），仅用于满足 GWASFlow BED 路线的身份对账契约；
- 效应方向以 GEMMA allele1(reported) 为准（GWASFlow 默认 target_effect_allele=reported），
  gwas.tsv 的 ref/alt 列语义为该登记约定，不构成参考等位基因声明。

支持两种输入（产物应逐字节一致，可交叉对账）：
  --tped <atpolydb.tped>   转换前（plink 未装时可用）
  --bim  <genotype.bim>    转换后（plink --keep-allele-order）

用法：
  python3 make_variant_map_from_bim.py --tped data/smoke/atpolydb.tped \
      --out data/smoke/variant_map.tsv
"""

import argparse
import csv
import sys
from pathlib import Path

VALID_BASES = {"A", "C", "G", "T"}
COLUMNS = ["variant_id", "original_id", "chrom", "engine_chrom", "pos", "ref", "alt"]


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="make_variant_map_from_bim.py",
        description="AtPolyDB TPED/BIM → GWASFlow variant_map.tsv（观测等位基因登记约定）",
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--tped", help="atpolydb.tped（转换前输入）")
    source.add_argument("--bim", help="genotype.bim（plink --keep-allele-order 转换后）")
    parser.add_argument("--chr-count", type=int, default=5,
                        help="染色体对数（拟南芥 5）")
    parser.add_argument("--out", required=True, help="variant_map.tsv 输出路径")
    return parser.parse_args(argv)


def die(message):
    print(f"[make_variant_map ERROR] {message}", file=sys.stderr)
    sys.exit(2)


def parse_tped(path):
    """TPED：chrom rsid cM pos + 每样本2列基因型（AtPolyDB 自交系纯合编码）。

    无显式等位基因列：从基因型列收集观测等位基因（"0" 为缺失）；
    恰好两种 → (ref, alt) = 字典序 (小, 大)（登记约定，集合与 BIM 对账一致）。
    """
    records = []
    seen = set()
    with open(path, encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6 or (len(parts) - 4) % 2 != 0:
                die(f"tped 第 {line_no} 行列结构异常（{len(parts)} 列）")
            chrom, rsid, _cm, pos = parts[:4]
            alleles = set()
            for i in range(4, len(parts), 2):
                for token in (parts[i], parts[i + 1]):
                    token = token.strip().upper()
                    if token and token != "0":
                        alleles.add(token)
            if len(alleles) != 2:
                continue   # 单态/超双等位：跳过（计数由调用方汇总）
            key = (chrom, pos)
            if key in seen:
                die(f"tped 第 {line_no} 行位点坐标重复: {key}")
            seen.add(key)
            ref, alt = sorted(alleles)   # 登记约定：ref=字典序较小
            records.append((chrom, rsid, pos, alt, ref))
    return records


def _count_tped_lines(path):
    with open(path, encoding="utf-8") as handle:
        return sum(1 for line in handle if line.strip())


def parse_bim(path):
    """BIM：chrom rsid cM pos a1 a2 → 同构记录（--keep-allele-order 下与 tped 一致）。"""
    records = []
    seen_ids = set()
    with open(path, encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            parts = line.split()
            if len(parts) < 6:
                die(f"bim 第 {line_no} 行列数不足（<6）")
            chrom, rsid, _cm, pos, a1, a2 = parts[:6]
            if rsid in seen_ids:
                die(f"bim 第 {line_no} 行 rsid 重复: {rsid}")
            seen_ids.add(rsid)
            records.append((chrom, rsid, pos, a1.upper(), a2.upper()))
    return records


def main(argv=None):
    args = parse_args(argv)

    if args.tped:
        records = parse_tped(args.tped)
        n_input_total = _count_tped_lines(args.tped)
        convention = "ref=字典序较小/alt=字典序较大（从基因型列推断观测等位基因）"
    else:
        records = parse_bim(args.bim)
        n_input_total = len(records)
        convention = "ref=A2/alt=A1（跟随 plink --keep-allele-order 列序）"

    if not records:
        die("输入无双等价可记录位点")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    n_bad = n_input_total - len(records)
    for chrom, rsid, pos, a1, a2 in records:
        # engine_chrom 校验（数值 1..chr_count；未知/越界即报错，不压 0）
        try:
            engine = int(chrom)
        except ValueError:
            die(f"染色体非数值（BED 路线要求 BIM 已是引擎数值染色体）: {chrom!r} @ {rsid}")
        if not (1 <= engine <= args.chr_count):
            die(f"染色体超出 1..{args.chr_count}: {engine} @ {rsid}")
        try:
            pos_int = int(pos)
        except ValueError:
            die(f"位置非整数: {pos!r} @ {rsid}")
        if pos_int <= 0:
            die(f"位置非正整数: {pos} @ {rsid}")
        if a1 not in VALID_BASES or a2 not in VALID_BASES or a1 == a2:
            n_bad += 1
            n_bad += 1
            continue   # 非双等位 A/C/G/T：记录并跳过（GWASFlow normalize BED 路线会对账拦截）
        # 登记约定见 convention（集合与 BIM 对账一致；方向为登记约定非参考断言）
        ref, alt = a2, a1
        variant_id = f"{chrom}:{pos}:{ref}:{alt}"
        rows.append([variant_id, rsid, str(chrom), str(engine), str(pos_int), ref, alt])

    # variant_id 唯一性
    ids = [r[0] for r in rows]
    if len(set(ids)) != len(ids):
        dupes = sorted({i for i in ids if ids.count(i) > 1})[:5]
        die(f"variant_id 重复（同坐标双等位冲突）: {dupes}")

    with open(out_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(COLUMNS)
        for row in rows:
            writer.writerow(row)

    print(
        f"[make_variant_map] {out_path}: {len(rows)} 位点"
        f"（输入 {n_input_total}，跳过 {n_bad}）；约定: {convention}；"
        f"观测等位基因，非参考方向断言"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
