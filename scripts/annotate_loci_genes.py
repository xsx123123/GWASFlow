#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GWASFlow 候选基因注释（参数手册 §6 契约重实现）。

GFF3 gene 区间相交（1-based inclusive、flank_bp 默认 0 且为已定义区域外的额外延伸）：
- 组装一致性强制（基因染色体与 loci 染色体完全无法对应 → 阻断，无"忽略匹配错误"开关）
- 功能表按唯一 gene_id 连接；一基因多行先汇总（禁止直接 join 造成重复）
- 无重叠基因保留 no_overlap 行（locus 级）
- 空 loci → 合法表头输出；旁写 .meta.json

错误码：0=成功；2=格式/配置/输入错误。
"""

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

EXIT_OK = 0
EXIT_INPUT_ERROR = 2

COLUMNS = [
    "locus_id", "chrom", "locus_start", "locus_end",
    "gene_id", "gene_start", "gene_end", "strand",
    "distance_to_locus_bp", "relation", "description", "status",
]


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="annotate_loci_genes.py",
        description="将已定义候选区域与 GFF3 gene 位置相交，按精确 gene_id 连接可选描述。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--loci", required=True, help="loci.tsv（locus_id,chrom,start,end，1-based inclusive）")
    parser.add_argument("--gff", required=True, help="未压缩 GFF3（默认读取 gene feature，不支持 GTF）")
    parser.add_argument("--genome-version", required=True, help="声明共同组装（不会自动证明序列一致）")
    parser.add_argument("--gene-feature", default="gene", help="GFF 第三列要匹配的特征类型")
    parser.add_argument("--id-attribute", default="ID", help="GFF3 属性中的唯一基因 ID")
    parser.add_argument("--flank-bp", type=int, default=0,
                        help="在输入区域外再扩展的 bp（>=0；默认 0 不扩大窗口）")
    parser.add_argument("--function-table", default=None,
                        help="功能表 TSV（gene_id 与 description 列）")
    parser.add_argument("--out", required=True, help="候选表输出路径（旁写 .meta.json）")
    return parser.parse_args(argv)


def fail(message, code=EXIT_INPUT_ERROR):
    print(f"[annotate_loci_genes ERROR] {message}", file=sys.stderr)
    sys.exit(code)


def sha256_of_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_attributes(text):
    attrs = {}
    for token in text.split(";"):
        if "=" in token:
            key, value = token.split("=", 1)
            attrs[key.strip()] = value.strip()
    return attrs


def read_gff_genes(gff_path, gene_feature, id_attribute):
    """读取 GFF3 gene feature：{chrom: [（gene_id, start, end, strand）]}。"""
    genes_by_chrom = {}
    n_features = 0
    with open(gff_path, encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != gene_feature:
                continue
            chrom = parts[0]
            try:
                start, end = int(parts[3]), int(parts[4])
            except ValueError:
                fail(f"GFF 第 {line_no} 行坐标非整数: {parts[3]}/{parts[4]}")
            if start > end:
                start, end = end, start
            strand = parts[6] if parts[6] in ("+", "-") else "."
            attrs = parse_attributes(parts[8])
            gene_id = attrs.get(id_attribute, "")
            if not gene_id:
                fail(
                    f"GFF 第 {line_no} 行 gene feature 缺失属性 '{id_attribute}'；"
                    f"请用 --id-attribute 指定正确属性"
                )
            genes_by_chrom.setdefault(chrom, []).append({
                "gene_id": gene_id, "start": start, "end": end, "strand": strand,
            })
            n_features += 1
    return genes_by_chrom, n_features


def read_function_table(path):
    """功能表：一基因多行先汇总（合并描述），禁止直接 join 造成候选行重复。"""
    if not path:
        return {}
    descriptions = {}
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = [c.strip() for c in (reader.fieldnames or [])]
        if "gene_id" not in fieldnames or "description" not in fieldnames:
            fail(f"功能表缺失列 gene_id/description（实际列: {fieldnames}）")
        for line_no, row in enumerate(reader, start=2):
            gene_id = (row.get("gene_id") or "").strip()
            desc = (row.get("description") or "").strip()
            if not gene_id:
                continue
            if gene_id in descriptions:
                descriptions[gene_id] = f"{descriptions[gene_id]}; {desc}"
            else:
                descriptions[gene_id] = desc
    return descriptions


def read_loci(path):
    loci = []
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            try:
                loci.append({
                    "locus_id": row["locus_id"],
                    "chrom": row["chrom"],
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                })
            except (KeyError, ValueError) as exc:
                fail(f"loci.tsv 行解析失败: {row}（{exc}）")
    return loci


def main(argv=None):
    args = parse_args(argv)

    if args.flank_bp < 0:
        fail(f"--flank-bp 必须 >=0（当前 {args.flank_bp}；默认 0 不扩大窗口）")

    loci = read_loci(args.loci)
    genes_by_chrom, n_features = read_gff_genes(args.gff, args.gene_feature, args.id_attribute)
    descriptions = read_function_table(args.function_table)

    # ---- 组装一致性强制（无"忽略匹配错误"开关）----
    if loci and n_features > 0:
        loci_chroms = {l["chrom"] for l in loci}
        gff_chroms = set(genes_by_chrom.keys())
        if not (loci_chroms & gff_chroms):
            fail(
                f"基因染色体与 loci 染色体完全无法对应（疑似命名/组装不一致）："
                f"loci 染色体样例 {sorted(loci_chroms)[:5]}，"
                f"GFF 染色体样例 {sorted(gff_chroms)[:5]}。"
                f"请检查 GFF3 的基因组版本/命名与当前组装（{args.genome_version}）是否一致；"
                f"当前适配器不提供'忽略匹配错误'开关"
            )

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for locus in loci:
        window_start = max(1, locus["start"] - args.flank_bp)
        window_end = locus["end"] + args.flank_bp
        genes = genes_by_chrom.get(locus["chrom"], [])
        in_window = [g for g in genes
                     if not (g["end"] < window_start or g["start"] > window_end)]
        if not in_window:
            rows.append([
                locus["locus_id"], locus["chrom"], locus["start"], locus["end"],
                "NA", "NA", "NA", "NA", "NA", "no_overlap", "NA", "no_overlap",
            ])
            continue
        for gene in sorted(in_window, key=lambda g: (g["start"], g["gene_id"])):
            if not (gene["end"] < locus["start"] or gene["start"] > locus["end"]):
                relation, distance = "overlap", 0
            elif gene["end"] < locus["start"]:
                relation = "upstream"
                distance = locus["start"] - gene["end"]
            else:
                relation = "downstream"
                distance = gene["start"] - locus["end"]
            rows.append([
                locus["locus_id"], locus["chrom"], locus["start"], locus["end"],
                gene["gene_id"], gene["start"], gene["end"], gene["strand"],
                distance, relation,
                descriptions.get(gene["gene_id"], "NA"), "ok",
            ])

    with open(out_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(COLUMNS)
        for row in rows:
            writer.writerow(row)

    # ---- .meta.json（来源元数据）----
    meta = {
        "schema": "gwasflow.candidates.v1",
        "genome_version": args.genome_version,
        "gene_feature": args.gene_feature,
        "id_attribute": args.id_attribute,
        "flank_bp": args.flank_bp,
        "coordinate_system": "1-based inclusive（distance 为基因与原 locus 最短坐标距离；重叠为 0）",
        "n_loci": len(loci),
        "n_gene_features": n_features,
        "n_candidate_rows": len(rows),
        "n_loci_no_overlap": sum(1 for r in rows if r[9] == "no_overlap"),
        "sources": {
            "loci": {"path": str(Path(args.loci).resolve()),
                     "sha256": sha256_of_file(args.loci)},
            "gff": {"path": str(Path(args.gff).resolve()),
                    "sha256": sha256_of_file(args.gff)},
        },
        "notes": (
            "不做 GO/KEGG 富集、转录本后果预测或因果推断；"
            "功能表一基因多行已按规则汇总（'; ' 连接）"
        ),
    }
    if args.function_table:
        meta["sources"]["function_table"] = {
            "path": str(Path(args.function_table).resolve()),
            "sha256": sha256_of_file(args.function_table),
        }
    meta_path = out_path.with_suffix(out_path.suffix + ".meta.json")
    with open(meta_path, "w", encoding="utf-8") as handle:
        json.dump(meta, handle, indent=2, ensure_ascii=False)

    print(
        f"[annotate_loci_genes] 候选注释完成：{len(rows)} 行"
        f"（{len(loci)} loci，flank={args.flank_bp}bp）"
    )
    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
