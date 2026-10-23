#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GWASFlow 显著区域定义（参数手册 §9.5）。

- lead 只在选定主模型与其阈值上选取（阈值从 association metadata 唯一读取，
  禁止另抄一份不同阈值）
- PLINK2 --clump LD 去冗余（clump 不等于精细定位）
- region-boundary=ld|fixed_window 显式标注（fixed_window 不称可信区间）
- locus 登记字段齐全（locus_id/trait/model/chrom/start/end/lead_variant/lead_p/
  boundary_method/参数/回退原因/ld_sample 摘要；start<=end、1-based inclusive）
- 区域定义使用该性状实际材料（keep）
- 无显著位点 → 合法空三表正常交付

错误码：0=成功；2=格式/配置/输入错误。
"""

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

EXIT_OK = 0
EXIT_INPUT_ERROR = 2


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="define_loci.py",
        description="主检验显著区域定义：lead 选取 + LD clump 去冗余 + 区域边界登记。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--sumstats", required=True, help="08.results/{trait}/gwas.tsv")
    parser.add_argument("--significant", required=True, help="08.results/{trait}/significant.tsv")
    parser.add_argument("--association-metadata", required=True,
                        help="08.results/{trait}/metadata.json（阈值唯一来源）")
    parser.add_argument("--trait", required=True, help="性状 ID")
    parser.add_argument("--model", default="gemma_lmm", help="主模型（登记用）")
    parser.add_argument("--bfile", required=True, help="cohort BED 前缀")
    parser.add_argument("--keep", required=True, help="该性状有效材料 keep.tsv")
    parser.add_argument("--mode", default="ld_clump", choices=["ld_clump"], help="去冗余模式")
    parser.add_argument("--p-threshold-from-metadata", action="store_true", default=True,
                        help="阈值从 metadata 读取（不再另抄一份不同阈值）")
    parser.add_argument("--clump-r2", type=float, default=0.20, help="clump r2")
    parser.add_argument("--clump-kb", type=int, default=500, help="clump 搜索半径（kb）")
    parser.add_argument("--region-boundary", default="ld", choices=["ld", "fixed_window"],
                        help="区域边界方法（fixed_window 仅标识方法，不称可信区间）")
    parser.add_argument("--fallback-flank-bp", type=int, default=50000,
                        help="LD 信息不足时回退窗口（bp，整数>=0）")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--out-dir", required=True, help="输出目录（新/空目录）")
    return parser.parse_args(argv)


def fail(message, code=EXIT_INPUT_ERROR):
    print(f"[define_loci ERROR] {message}", file=sys.stderr)
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


def read_gwas_rows(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def parse_decimal(text):
    from decimal import Decimal, InvalidOperation
    try:
        return Decimal(str(text).strip())
    except (InvalidOperation, ValueError):
        return None


def parse_clumped(path):
    """解析 plink2 .clumped：返回 {index_snp: [member_snps]}。"""
    if not Path(path).exists():
        return None
    clusters = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            index_snp = parts[1] if not parts[0].isdigit() else parts[1]
            # plink2 .clumped: CHROM SNP ... SP2(成员=值分号分隔)
            members = []
            if len(parts) >= 10:
                sp2 = parts[9] if "SP2" in line.upper() else ""
                for token in sp2.split(";"):
                    if token and not token.startswith("NONE"):
                        member = token.split("=")[0].strip()
                        if member:
                            members.append(member)
            clusters[index_snp] = members
    return clusters


def main(argv=None):
    args = parse_args(argv)
    out_dir = ensure_empty_dir(args.out_dir)

    # ---- 阈值唯一来源：association metadata ----
    with open(args.association_metadata, encoding="utf-8") as handle:
        assoc_meta = json.load(handle)
    mt = assoc_meta.get("multiple_testing", {})
    threshold_text = mt.get("threshold", "NA")

    if assoc_meta.get("status") == "failed_software":
        _write_empty_loci(args, out_dir, assoc_meta,
                          "上游关联分析失败（failed_software）；无显著区域")
        print("[define_loci] 上游失败传导：输出合法空三表")
        return EXIT_OK

    threshold = parse_decimal(threshold_text)
    if threshold is None:
        fail(
            f"association metadata 阈值不可用: {threshold_text!r}"
            f"（状态 {assoc_meta.get('status')!r}）"
        )

    # ---- 显著位点读取 ----
    sig_rows = read_gwas_rows(args.significant)
    if not sig_rows:
        _write_empty_loci(args, out_dir, assoc_meta, "无显著位点（合法空结果）")
        print("[define_loci] 无显著位点：合法空三表正常交付")
        return EXIT_OK

    # ---- clump 输入表（variant_id + p_value，含 chrom/pos 辅助列）----
    clump_input = out_dir / "clump_input.tsv"
    with open(clump_input, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["variant_id", "p_value"])
        for row in sig_rows:
            p_text = row.get("p_value", "")
            p = parse_decimal(p_text)
            if p is not None and p != 0:
                writer.writerow([row["variant_id"], format(p, "g")])
            elif row.get("status") == "p_underflow":
                # 下溢 P 以可表示极小值参与 clump 排序（不改统计表原值）
                writer.writerow([row["variant_id"], "1e-300"])

    # ---- PLINK2 clump 去冗余（该性状实际材料）----
    log_path = out_dir / "clump.log"
    clump_prefix = out_dir / "clumped"
    argv_cmd = [
        "plink2", "--bfile", args.bfile,
        "--keep", args.keep,
        "--clump", str(clump_input),
        "--clump-snp-field", "variant_id",
        "--clump-p-field", "p_value",
        "--clump-r2", str(args.clump_r2),
        "--clump-kb", str(args.clump_kb),
        "--threads", str(args.threads),
        "--out", str(clump_prefix),
    ]
    clump_ok = True
    with open(log_path, "w", encoding="utf-8") as handle:
        handle.write(f"=== RUN: {' '.join(argv_cmd)} ===\n")
        handle.flush()
        result = subprocess.run(argv_cmd, stdout=handle, stderr=subprocess.STDOUT,
                                check=False)
    if result.returncode != 0:
        clump_ok = False

    clusters = parse_clumped(f"{clump_prefix}.clumped") if clump_ok else None
    fallback_reason = ""
    if not clusters:
        fallback_reason = (
            "clump 去冗余不可用（plink2 退出码 "
            f"{result.returncode if not clump_ok else 'parse'}）；全部显著位点保留为 lead"
        )

    # ---- lead 集合与区域定义 ----
    keep_sha = sha256_of_file(args.keep)
    loci_rows = []
    lead_rows = []

    def chrom_sort_key(chrom):
        try:
            return int(chrom)
        except ValueError:
            return 10**9

    sig_sorted = sorted(sig_rows, key=lambda r: (chrom_sort_key(r.get("chrom", "")),
                                                 int(r.get("pos", 0))))
    used_leads = set(clusters.keys()) if clusters else set()

    for row in sig_sorted:
        vid = row["variant_id"]
        is_lead = vid in used_leads if clusters else True
        if not is_lead:
            continue
        chrom = row.get("chrom", "")
        try:
            pos = int(row.get("pos", 0))
        except ValueError:
            continue
        p_text = row.get("p_value", "")

        members = clusters.get(vid, []) if clusters else []
        boundary_method = args.region_boundary
        params = json.dumps({"clump_r2": args.clump_r2, "clump_kb": args.clump_kb})
        fb_reason = fallback_reason

        if args.region_boundary == "ld" and members:
            member_pos = []
            for member_vid in members:
                for sig_row in sig_sorted:
                    if sig_row["variant_id"] == member_vid:
                        try:
                            member_pos.append(int(sig_row.get("pos", 0)))
                        except ValueError:
                            pass
                        break
            if member_pos:
                start, end = min(pos, *member_pos), max(pos, *member_pos)
            else:
                start, end = max(1, pos - args.fallback_flank_bp), pos + args.fallback_flank_bp
                boundary_method = "fixed_window"
                fb_reason = fb_reason or "LD 成员坐标不可解析，回退固定窗口（仅方法标识，非可信区间）"
        elif args.region_boundary == "ld":
            start = max(1, pos - args.fallback_flank_bp)
            end = pos + args.fallback_flank_bp
            if not fb_reason:
                fb_reason = "无 LD 成员信息，回退固定窗口（仅方法标识，非可信区间）"
        else:
            start = max(1, pos - args.fallback_flank_bp)
            end = pos + args.fallback_flank_bp

        locus_id = f"locus_{len(loci_rows) + 1:04d}"
        lead_rows.append({
            "locus_id": locus_id, "variant_id": vid, "chrom": chrom, "pos": pos,
            "p_value": p_text, "p_status": row.get("status", ""),
        })
        loci_rows.append({
            "locus_id": locus_id, "trait": args.trait, "model": args.model,
            "chrom": chrom, "start": start, "end": end,
            "lead_variant": vid, "lead_p": p_text,
            "boundary_method": boundary_method, "boundary_params": params,
            "fallback_reason": fb_reason,
            "ld_sample_sha256": keep_sha,
            "n_clump_members": len(members),
        })

    _write_loci_outputs(args, out_dir, assoc_meta, lead_rows, loci_rows,
                        threshold_text, keep_sha)
    print(
        f"[define_loci] 区域定义完成：{len(lead_rows)} 个 lead/locus"
        f"（阈值 {threshold_text}，boundary={args.region_boundary}）"
    )
    return EXIT_OK


def _write_loci_outputs(args, out_dir, assoc_meta, lead_rows, loci_rows,
                        threshold_text, keep_sha):
    with open(out_dir / "lead.tsv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["locus_id", "variant_id", "chrom", "pos", "p_value", "p_status"])
        for r in lead_rows:
            writer.writerow([r["locus_id"], r["variant_id"], r["chrom"],
                             r["pos"], r["p_value"], r["p_status"]])
    with open(out_dir / "loci.tsv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "locus_id", "trait", "model", "chrom", "start", "end",
            "lead_variant", "lead_p", "boundary_method", "boundary_params",
            "fallback_reason", "ld_sample_sha256", "n_clump_members",
        ])
        for r in loci_rows:
            writer.writerow([
                r["locus_id"], r["trait"], r["model"], r["chrom"], r["start"], r["end"],
                r["lead_variant"], r["lead_p"], r["boundary_method"],
                r["boundary_params"], r["fallback_reason"], r["ld_sample_sha256"],
                r["n_clump_members"],
            ])
    manifest = {
        "schema": "gwasflow.loci.v1",
        "trait": args.trait,
        "model": args.model,
        "p_threshold": threshold_text,
        "p_threshold_source": "association metadata（唯一来源，未另抄）",
        "clump_params": {"r2": args.clump_r2, "kb": args.clump_kb, "mode": args.mode},
        "region_boundary": args.region_boundary,
        "fallback_flank_bp": args.fallback_flank_bp,
        "ld_sample_sha256": keep_sha,
        "n_loci": len(loci_rows),
        "note": "clump 用于去冗余而非精细定位；fixed_window 仅方法标识，不称可信区间",
    }
    with open(out_dir / "locus_manifest.json", "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)


def _write_empty_loci(args, out_dir, assoc_meta, reason):
    _write_loci_outputs(args, out_dir, assoc_meta, [], [],
                        str(assoc_meta.get("multiple_testing", {}).get("threshold", "NA")),
                        sha256_of_file(args.keep))
    manifest_path = out_dir / "locus_manifest.json"
    with open(manifest_path, encoding="utf-8") as handle:
        manifest = json.load(handle)
    manifest["empty_reason"] = reason
    with open(manifest_path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)


if __name__ == "__main__":
    sys.exit(main())
