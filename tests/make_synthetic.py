#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GWASFlow 合成测试数据生成器（L1/L2/L3 测试共用）。

生成自洽小群体数据集：
- samples.csv（n 材料）、phenotypes.tsv（两性状不同缺失模式）、traits.tsv、covariates.tsv
- BED/BIM/FAM 三件套（硬分型双等位 SNP，9 染色体 x snps_per_chrom 位点）
- VCF.gz + 索引占位（供门控组件校验；bcftools 真实转换由 L3 mock 驱动）
- 染色体重命名映射/染色体表、GFF3（含 no_overlap 基因）、功能表（一基因多行）
- FASTA + FAI（REF 校验基准）
- GEMMA 风格 assoc 表、PCA eigenvec、K 矩阵、run_metadata（表级注入测试用）

边界数据构造器（--boundary）：
  multiallelic / ref_mismatch / unknown_chrom / coord_mismatch /
  p_underflow / collinear / zero_significant / duplicate_fid_iid / missing_index

用法:
  python3 tests/make_synthetic.py --out-dir /tmp/gwas_data [--n-samples 30]
      [--snps-per-chrom 40] [--boundary NAME]
"""

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import random
import sys
from pathlib import Path

CHROMS = [f"Lsat{i}" for i in range(1, 10)]   # 生菜 9 染色体占位（M2）
ENGINE_CHROMS = {c: i + 1 for i, c in enumerate(CHROMS)}


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="GWASFlow 合成测试数据生成器")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--n-samples", type=int, default=30)
    parser.add_argument("--snps-per-chrom", type=int, default=40)
    parser.add_argument("--seed", type=int, default=20261023)
    parser.add_argument("--boundary", default=None,
                        choices=[None, "multiallelic", "ref_mismatch", "unknown_chrom",
                                 "coord_mismatch", "p_underflow", "collinear",
                                 "zero_significant", "duplicate_fid_iid", "missing_index"],
                        help="边界数据构造器")
    return parser.parse_args(argv)


def write_text(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)


def md5_of_file(path):
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        digest.update(handle.read())
    return digest.hexdigest()


def generate_variants(args, rng):
    """合成双等位 SNP 列表：[(chrom, pos, ref, alt)]。"""
    variants = []
    bases = ["A", "C", "G", "T"]
    for chrom in CHROMS:
        positions = sorted(rng.sample(range(1000, 900000), args.snps_per_chrom))
        for pos in positions:
            ref = rng.choice(bases)
            alt = rng.choice([b for b in bases if b != ref])
            variants.append((chrom, pos, ref, alt))
    return variants


def genotype_matrix(args, rng, variants):
    """硬分型基因型矩阵：{variant_idx: [gt per sample]}（0/1/2/None=缺失）。"""
    matrix = []
    for _ in variants:
        maf = rng.uniform(0.05, 0.5)
        row = []
        for _ in range(args.n_samples):
            if rng.random() < 0.02:   # 2% 缺失
                row.append(None)
            else:
                hom = rng.random() < maf * maf
                het = rng.random() < 2 * maf * (1 - maf)
                row.append(1 if het else (2 if hom else 0))
        matrix.append(row)
    return matrix


def write_bed(prefix, samples, variants, matrix):
    """写 PLINK BED/BIM/FAM 三件套（SNP-major）。"""
    n = len(samples)
    with open(f"{prefix}.fam", "w", encoding="utf-8") as handle:
        for s in samples:
            handle.write(f"{s}\t{s}\t0\t0\t0\t-9\n")
    with open(f"{prefix}.bim", "w", encoding="utf-8") as handle:
        for i, (chrom, pos, ref, alt) in enumerate(variants):
            vid = f"{chrom}:{pos}:{ref}:{alt}"
            handle.write(f"{ENGINE_CHROMS[chrom]}\t{vid}\t0\t{pos}\t{alt}\t{ref}\n")
    with open(f"{prefix}.bed", "wb") as handle:
        handle.write(bytes([0x6C, 0x1B, 0x01]))
        for row in matrix:
            nbytes = (n + 3) // 4
            buf = bytearray(nbytes)
            for j, gt in enumerate(row):
                if gt is None:
                    code = 0
                elif gt == 1:
                    code = 1
                elif gt == 2:
                    code = 3
                else:
                    code = 2
                byte_idx, bit_idx = j // 4, (j % 4) * 2
                buf[byte_idx] |= code << bit_idx
            handle.write(bytes(buf))


def write_vcf(path, samples, variants, matrix, boundary=None, ref_override=None):
    """写 VCF.gz（硬分型 GT）+ 触碰 .csi 索引占位。"""
    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=Lsat1,length=1000000>",
    ]
    lines += [f"##contig=<ID={c},length=1000000>" for c in CHROMS]
    lines += ['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)]
    extra = []
    gt_all_ref = "\t".join("0/0" for _ in samples)
    if boundary == "multiallelic":
        extra.append(f"Lsat1\t999999\t.\tA\tG,T\t25\tPASS\t.\tGT\t{gt_all_ref}")
        extra.append(f"Lsat1\t999998\t.\tAT\tGC\t25\tPASS\t.\tGT\t{gt_all_ref}")
    if boundary == "unknown_chrom":
        extra.append(f"ChrUnplaced\t500\t.\tA\tG\t25\tPASS\t.\tGT\t{gt_all_ref}")
    if boundary == "ref_mismatch":
        # REF 与 FASTA 不一致（FASTA 中该位置是 G，VCF 写 A）
        extra.append(f"Lsat2\t555001\t.\tA\tC\t25\tPASS\t.\tGT\t{gt_all_ref}")

    for i, (chrom, pos, ref, alt) in enumerate(variants):
        ref_here = ref_override((chrom, pos)) if ref_override and ref_override((chrom, pos)) else ref
        vid = "."
        gts = []
        for gt in matrix[i]:
            gts.append("./." if gt is None else {0: "0/0", 1: "0/1", 2: "1/1"}[gt])
        lines.append(f"{chrom}\t{pos}\t{vid}\t{ref_here}\t{alt}\t25\tPASS\t.\tGT\t" + "\t".join(gts))
    lines.extend(extra)

    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    Path(str(path) + ".csi").touch()   # 索引占位（门控组件校验用）


def write_fasta(path, variants):
    """参考 FASTA：每染色体固定种子序列 + 变异位点写入 REF 碱基。"""
    per_chrom = {c: {} for c in CHROMS}
    for chrom, pos, ref, alt in variants:
        per_chrom[chrom][pos] = ref
    fai_lines = []
    with open(path, "w", encoding="utf-8") as handle:
        for chrom in CHROMS:
            length = 1000000
            seq = ["A"] * length
            for pos, base in per_chrom[chrom].items():
                if pos <= length:
                    seq[pos - 1] = base
            handle.write(f">{chrom}\n")
            for i in range(0, length, 60):
                handle.write("".join(seq[i:i + 60]) + "\n")
            fai_lines.append(f"{chrom}\t{length}\t{offset_of(path, chrom)}\t60\t61\n")
    with open(f"{path}.fai", "w", encoding="utf-8") as handle:
        handle.writelines(fai_lines)


def offset_of(path, chrom):
    return 0   # FAI offset 仅占位（bcftools mock 场景不真实读取）


def main(argv=None):
    args = parse_args(argv)
    rng = random.Random(args.seed)
    out = Path(args.out_dir).resolve()
    if out.exists() and any(out.iterdir()):
        print(f"[make_synthetic ERROR] 输出目录非空: {out}", file=sys.stderr)
        sys.exit(2)
    out.mkdir(parents=True)

    n = args.n_samples
    samples = [f"S{i:03d}" for i in range(1, n + 1)]
    variants = generate_variants(args, rng)
    matrix = genotype_matrix(args, rng, variants)

    # ---- samples.csv ----
    if args.boundary == "duplicate_fid_iid":
        rows = ["sample,fid,iid,sample_name,group,include"]
        for i, s in enumerate(samples):
            rows.append(f"{s},{s},{s},name_{s},G{i % 3},true")
        rows.append(f"{samples[0]}_dup,{samples[0]},{samples[0]},dup,G0,true")
        write_text(out / "samples.csv", "\n".join(rows) + "\n")
    else:
        rows = ["sample,fid,iid,sample_name,group,include"]
        for i, s in enumerate(samples):
            rows.append(f"{s},{s},{s},name_{s},G{i % 3},true")
        write_text(out / "samples.csv", "\n".join(rows) + "\n")

    # ---- phenotypes.tsv（两性状不同缺失模式：trait1 缺 S005/S011/S017；trait2 缺 S003/S008/S013/S019/S024/S029）----
    missing_t1 = {"S005", "S011", "S017"}
    missing_t2 = {"S003", "S008", "S013", "S019", "S024", "S029"}
    rows = ["sample\tleaf_length\tleaf_width"]
    for i, s in enumerate(samples):
        v1 = "NA" if s in missing_t1 else f"{rng.uniform(5, 25):.2f}"
        v2 = "NA" if s in missing_t2 else f"{rng.uniform(2, 12):.2f}"
        if args.boundary == "zero_significant" or True:
            # 表型与基因型独立（默认无真实关联 → 零显著路径可复现）
            pass
        rows.append(f"{s}\t{v1}\t{v2}")
    write_text(out / "phenotypes.tsv", "\n".join(rows) + "\n")

    # ---- traits.tsv ----
    covar_num = "" if args.boundary == "collinear" else "days_to_measure"
    write_text(out / "traits.tsv",
               "trait\tphenotype_column\ttrait_type\tunit\ttransform\tprimary_model"
               "\tnumeric_covariates\tcategorical_covariates\tn_pcs\n"
               f"leaf_length\tleaf_length\tcontinuous\tcm\tnone\tgemma_lmm\t{covar_num}\tbatch\t3\n"
               "leaf_width\tleaf_width\tcontinuous\tcm\tnone\tgemma_lmm\t\t\t0\n")

    # ---- covariates.tsv ----
    if args.boundary == "collinear":
        rows = ["sample\tdays_to_measure\tdays_copy\tbatch"]
        for i, s in enumerate(samples):
            d = rng.uniform(10, 30)
            rows.append(f"{s}\t{d:.2f}\t{d:.2f}\tB{1 + i % 2}")
        write_text(out / "covariates.tsv", "\n".join(rows) + "\n")
    else:
        rows = ["sample\tdays_to_measure\tbatch"]
        for i, s in enumerate(samples):
            rows.append(f"{s}\t{rng.uniform(10, 30):.2f}\tB{1 + i % 2}")
        write_text(out / "covariates.tsv", "\n".join(rows) + "\n")

    # ---- 染色体映射/表 ----
    write_text(out / "chromosome_rename.tsv",
               "\n".join(f"{c}\t{ENGINE_CHROMS[c]}" for c in CHROMS) + "\n")
    write_text(out / "chromosome_table.tsv",
               "chrom\tlength_bp\n" + "\n".join(f"{c}\t1000000" for c in CHROMS) + "\n")

    # ---- 基因型三件套 + VCF ----
    write_bed(str(out / "cohort"), samples, variants, matrix)
    if args.boundary == "missing_index":
        write_vcf(out / "cohort.vcf.gz", samples, variants, matrix)
        (out / "cohort.vcf.gz.csi").unlink(missing_ok=True)
    else:
        write_vcf(out / "cohort.vcf.gz", samples, variants, matrix, boundary=args.boundary)

    # ---- FASTA + FAI（REF 校验基准）----
    write_fasta(out / "genome.fasta", variants)

    # ---- GFF3（gene feature；含与 loci 无关染色体基因 → no_overlap 路径）----
    gff_lines = ["##gff-version 3"]
    lead_chrom, lead_pos = variants[len(variants) // 2][0], variants[len(variants) // 2][1]
    for k in range(1, 21):
        chrom = CHROMS[k % len(CHROMS)]
        start = 100000 + k * 3000
        gff_lines.append(
            f"{chrom}\ttest\tgene\t{start}\t{start + 1500}\t.\t+\t.\tID=gene_{k:03d}"
        )
    write_text(out / "genes.gff3", "\n".join(gff_lines) + "\n")

    # ---- 功能表（gene_001 一基因多行 → 汇总去重演示）----
    write_text(out / "functions.tsv",
               "gene_id\tdescription\n"
               "gene_001\tresponse to auxin\n"
               "gene_001\tcell wall organization\n"
               "gene_002\tchloroplast development\n")

    # ---- PCA eigenvec（plink2 scols=fid 格式）----
    pca_lines = ["#FID\tIID\tPC1\tPC2\tPC3"]
    for i, s in enumerate(samples):
        pca_lines.append(
            f"{s}\t{s}\t{rng.gauss(0, 1):.4f}\t{rng.gauss(0, 1):.4f}\t{rng.gauss(0, 1):.4f}"
        )
    write_text(out / "pca.eigenvec", "\n".join(pca_lines) + "\n")

    # ---- GEMMA 风格 assoc 表（表级注入测试）----
    assoc_lines = ["chr\trs\tn_miss\tn_obs\tallele1\tallele2\taf\tbeta\tse\tlogl_H1\tl_remle\tp_wald"]
    for i, (chrom, pos, ref, alt) in enumerate(variants):
        vid = f"{chrom}:{pos}:{ref}:{alt}"
        n_miss = sum(1 for g in matrix[i] if g is None)
        n_obs = n - n_miss
        af = rng.uniform(0.05, 0.5)
        beta = rng.gauss(0, 0.3)
        se = abs(rng.gauss(0.05, 0.02))
        # 模拟 GEMMA allele1/allele2 方向不定（约半数行 allele1=ref）
        if i % 2 == 0:
            a1, a2 = alt, ref
        else:
            a1, a2 = ref, alt
        if args.boundary == "p_underflow" and i == 0:
            p_val = "0"
        elif args.boundary == "p_underflow" and i == 1:
            p_val = "1e-1000"
        elif args.boundary == "coord_mismatch" and i == 2:
            # 坐标错配（chrom 写错）
            assoc_lines.append(
                f"{ENGINE_CHROMS[chrom] + 1}\t{vid}\t{n_miss}\t{n_obs}\t{a1}\t{a2}"
                f"\t{af:.4f}\t{beta:.4f}\t{se:.4f}\t-100\t-100\t0.5"
            )
            continue
        elif args.boundary == "zero_significant" or True:
            p_val = f"{rng.uniform(0.05, 1):.6g}"
        assoc_lines.append(
            f"{ENGINE_CHROMS[chrom]}\t{vid}\t{n_miss}\t{n_obs}\t{a1}\t{a2}"
            f"\t{af:.4f}\t{beta:.4f}\t{se:.4f}\t-100\t-100\t{p_val}"
        )
    write_text(out / "association.assoc.txt", "\n".join(assoc_lines) + "\n")

    # ---- trait_mac.tsv（精确 MAC 注入测试）----
    mac_lines = ["variant_id\tchrom\tpos\tref\talt\talt_count\tmac\tallowed"]
    for i, (chrom, pos, ref, alt) in enumerate(variants):
        vid = f"{chrom}:{pos}:{ref}:{alt}"
        alt_count = sum(1 for g in matrix[i] if g == 1) + 2 * sum(1 for g in matrix[i] if g == 2)
        mac = min(alt_count, 2 * n - alt_count)
        mac_lines.append(
            f"{vid}\t{chrom}\t{pos}\t{ref}\t{alt}\t{alt_count}\t{mac}\t"
            f"{'true' if mac >= 5 else 'false'}"
        )
    write_text(out / "trait_mac.tsv", "\n".join(mac_lines) + "\n")

    # ---- run_metadata / sample metadata（harmonize 前置校验用）----
    sample_meta = {
        "schema": "gwasflow.trait_prep.v1", "trait": "leaf_length",
        "phenotype_column": "leaf_length", "transform": "none", "unit": "cm",
        "effect_scale": "cm", "n_total": n, "n_valid": n - 3,
        "primary_model": "gemma_lmm",
    }
    write_text(out / "sample_metadata.json", json.dumps(sample_meta, indent=2))

    run_meta = {
        "schema": "gwasflow.gemma_run.v1", "trait": "leaf_length", "test": "wald",
        "lmm_code": "1", "status": "software_exited_zero", "exit_code": 0,
        "bindings": {
            "trait_metadata": {
                "path": str(out / "sample_metadata.json"),
                "sha256": hashlib.sha256(
                    json.dumps(sample_meta, indent=2).encode()).hexdigest(),
                "n_valid": n - 3, "n_total": n, "trait": "leaf_length",
            },
        },
        "software_filters": {"maf": 0.05, "miss": 0.10, "hwe": 0, "r2": 0.9999},
    }
    write_text(out / "run_metadata.json", json.dumps(run_meta, indent=2))

    # ---- variant_map（harmonize 身份对账基准）----
    vm_lines = ["variant_id\toriginal_id\tchrom\tengine_chrom\tpos\tref\talt"]
    for chrom, pos, ref, alt in variants:
        vid = f"{chrom}:{pos}:{ref}:{alt}"
        vm_lines.append(f"{vid}\t.\t{chrom}\t{ENGINE_CHROMS[chrom]}\t{pos}\t{ref}\t{alt}")
    write_text(out / "variant_map.tsv", "\n".join(vm_lines) + "\n")

    # ---- K 矩阵 + sample_order（validate_kinship 测试）----
    k_matrix = [[1.0 if i == j else rng.uniform(-0.05, 0.15) for j in range(n)]
                for i in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            k_matrix[j][i] = k_matrix[i][j]   # 精确对称
    write_text(out / "kinship.cXX.txt",
               "\n".join("\t".join(f"{v:.6f}" for v in row) for row in k_matrix) + "\n")
    write_text(out / "sample_order.tsv",
               "fid\tiid\n" + "\n".join(f"{s}\t{s}" for s in samples) + "\n")

    # ---- expected_md5.tsv（提供方基线）----
    write_text(out / "expected_md5.tsv",
               "file\tmd5\n" +
               f"cohort.vcf.gz\t{md5_of_file(out / 'cohort.vcf.gz')}\n")

    # ---- 清单 ----
    print(f"[make_synthetic] 数据集生成完成: {out}")
    print(f"  材料 {n}，染色体 {len(CHROMS)}，位点 {len(variants)}")
    if args.boundary:
        print(f"  边界构造: {args.boundary}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
