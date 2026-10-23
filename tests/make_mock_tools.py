#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GWASFlow L3 测试用 mock 外部软件生成器。

生成假 plink2/gemma/bcftools/Rscript 可执行（输u5728预置合规文件），
驱动真实规则链端到端运行。mock 仅用于行为验证，不替代真实软件验收（M3）。
"""

import argparse
import stat
import sys
from pathlib import Path

PLINK2_MOCK = r'''#!/usr/bin/env python3
# mock plink2（GWASFlow L3 行为测试）：识别子命令并输出合规预置文件
import gzip
import random
import sys
from pathlib import Path

args = sys.argv[1:]
rng = random.Random(20261023)


def find_input():
    for i, a in enumerate(args):
        if a == "--bfile" and i + 1 < len(args):
            return ("bfile", args[i + 1])
    for i, a in enumerate(args):
        if a == "--vcf" and i + 1 < len(args):
            return ("vcf", args[i + 1])
    return (None, None)


input_kind, input_path = find_input()


def read_fam(prefix=None):
    rows = []
    if input_kind == "vcf":
        opener = gzip.open if str(input_path).endswith(".gz") else open
        with opener(input_path, "rt", errors="replace") as h:
            for line in h:
                if line.startswith("#CHROM"):
                    for s in line.strip().split("\t")[9:]:
                        rows.append([s, s, "0", "0", "0", "-9"])
                    break
    else:
        p = Path((prefix or input_path) + ".fam")
        if p.exists():
            for line in p.read_text().splitlines():
                parts = line.split()
                if len(parts) >= 2:
                    rows.append(parts)
    return rows


def read_bim(prefix=None):
    rows = []
    if input_kind == "vcf":
        opener = gzip.open if str(input_path).endswith(".gz") else open
        with opener(input_path, "rt", errors="replace") as h:
            for line in h:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 5:
                    chrom, pos, vid, ref, alt = parts[:5]
                    alts = alt.split(",")
                    stable_id = vid if vid not in (".", "", None) else \
                        f"{chrom}:{pos}:{ref}:{alts[0]}"
                    rows.append([chrom, stable_id, "0", pos, alts[0], ref])
    else:
        p = Path((prefix or input_path) + ".bim")
        if p.exists():
            for line in p.read_text().splitlines():
                parts = line.split()
                if len(parts) >= 6:
                    rows.append(parts)
    return rows


out_prefix = None
for i, a in enumerate(args):
    if a == "--out" and i + 1 < len(args):
        out_prefix = args[i + 1]

if out_prefix:
    fam = read_fam() if input_path else []
    bim = read_bim() if input_path else []

    if "--make-bed" in args:
        with open(out_prefix + ".fam", "w") as h:
            for r in fam:
                h.write("\t".join(r) + "\n")
        with open(out_prefix + ".bim", "w") as h:
            for r in bim:
                h.write("\t".join(r) + "\n")
        n = len(fam)
        with open(out_prefix + ".bed", "wb") as h:
            h.write(bytes([0x6C, 0x1B, 0x01]))
            for _ in bim:
                h.write(bytes((n + 3) // 4))
        with open(out_prefix + ".log", "w") as h:
            h.write("mock plink2 --make-bed done\n")
            h.write("3 variants removed due to missing genotype frequency (mock)\n")
    elif "--freq" in args:
        keep_n = len(fam)
        if "--keep" in args:
            ki = args.index("--keep")
            kf = args[ki + 1]
            if Path(kf).exists():
                keep_n = sum(1 for l in Path(kf).read_text().splitlines() if l.strip())
        with open(out_prefix + ".afreq", "w") as h:
            h.write("#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n")
            for r in bim:
                chrom, vid, _cm, pos, a1, a2 = r[0], r[1], r[2], r[3], r[4], r[5]
                af = rng.uniform(0.05, 0.5)
                h.write(f"{chrom}\t{vid}\t{a2}\t{a1}\t{af:.4f}\t{2 * keep_n}\n")
        with open(out_prefix + ".log", "w") as h:
            h.write("mock plink2 --freq done\n")
    elif "--missing" in args:
        with open(out_prefix + ".imiss", "w") as h:
            h.write("#IID\tMISS_PHENO\tN_MISS\tN_GENO\tF_MISS\n")
            for r in fam:
                h.write(f"{r[1]}\tY\t0\t{len(bim)}\t0.0\n")
        with open(out_prefix + ".lmiss", "w") as h:
            h.write("#CHROM\tID\tN_MISS\tN_GENO\tF_MISS\n")
            for r in bim:
                h.write(f"{r[0]}\t{r[1]}\t0\t{len(fam)}\t0.0\n")
        with open(out_prefix + ".log", "w") as h:
            h.write("mock plink2 --missing done\n")
    elif "--indep-pairwise" in args:
        keep = bim[: max(1, len(bim) * 2 // 3)]
        drop = bim[max(1, len(bim) * 2 // 3):]
        with open(out_prefix + ".prune.in", "w") as h:
            for r in keep:
                h.write(r[1] + "\n")
        with open(out_prefix + ".prune.out", "w") as h:
            for r in drop:
                h.write(r[1] + "\n")
        with open(out_prefix + ".log", "w") as h:
            h.write("mock plink2 --indep-pairwise done\n")
    elif "--pca" in args:
        n_pcs = 10
        for i, a in enumerate(args):
            if a == "--pca" and i + 1 < len(args) and args[i + 1].isdigit():
                n_pcs = int(args[i + 1])
        with open(out_prefix + ".eigenvec", "w") as h:
            h.write("#FID\tIID\t" + "\t".join(f"PC{k+1}" for k in range(n_pcs)) + "\n")
            for r in fam:
                vals = "\t".join(f"{rng.gauss(0, 1):.4f}" for _ in range(n_pcs))
                h.write(f"{r[0]}\t{r[1]}\t{vals}\n")
        with open(out_prefix + ".eigenval", "w") as h:
            for k in range(n_pcs):
                h.write(f"{rng.uniform(0.5, 3):.4f}\n")
        with open(out_prefix + ".log", "w") as h:
            h.write("mock plink2 --pca done\n")
    elif "--clump" in args:
        clump_file = None
        for i, a in enumerate(args):
            if a == "--clump" and i + 1 < len(args):
                clump_file = args[i + 1]
        rows = []
        if clump_file and Path(clump_file).exists():
            lines = Path(clump_file).read_text().splitlines()[1:]
            rows = [l.split("\t") for l in lines if l.strip()]
        with open(out_prefix + ".clumped", "w") as h:
            h.write("#CHROM\tSNP\tBP\tP\tTOTAL\tSP2\n")
            for r in rows[: max(1, len(rows) // 2)]:
                members = ";".join(f"{x[0]}={x[1]}" for x in rows[1:3])
                h.write(f"1\t{r[0]}\t100\t{r[1]}\t2\t{members or 'NONE'}\n")
        with open(out_prefix + ".log", "w") as h:
            h.write("mock plink2 --clump done\n")
    elif "--r2" in args:
        with open(out_prefix + ".ld", "w") as h:
            h.write("CHR_A\tSNP_A\tCHR_B\tSNP_B\tR2\n")
            if bim:
                for i in range(min(5, len(bim))):
                    a, b = bim[i % len(bim)], bim[(i + 1) % len(bim)]
                    h.write(f"{a[0]}\t{a[1]}\t{b[0]}\t{b[1]}\t0.8\n")
        with open(out_prefix + ".log", "w") as h:
            h.write("mock plink2 --r2 done\n")
    else:
        with open(out_prefix + ".log", "w") as h:
            h.write("mock plink2 (no-op)\n")

sys.exit(0)
''' 

GEMMA_MOCK = r'''#!/usr/bin/env python3
# mock gemma（GWASFlow L3 行为测试）
import os
import random
import sys
from pathlib import Path

args = sys.argv[1:]
rng = random.Random(20261023)

outdir = "."
outname = "output"
bfile = None
for i, a in enumerate(args):
    if a == "-outdir" and i + 1 < len(args):
        outdir = args[i + 1]
    if a == "-o" and i + 1 < len(args):
        outname = args[i + 1]
    if a == "-bfile" and i + 1 < len(args):
        bfile = args[i + 1]

if not args:
    print("gemma 0.98.1 (mock for GWASFlow L3 tests)\nUsage: gemma [options]")
    sys.exit(0)

def read_fam(prefix):
    rows = []
    p = Path(prefix + ".fam")
    if p.exists():
        for line in p.read_text().splitlines():
            parts = line.split()
            if len(parts) >= 2:
                rows.append((parts[0], parts[1]))
    return rows

def read_bim(prefix):
    rows = []
    p = Path(prefix + ".bim")
    if p.exists():
        for line in p.read_text().splitlines():
            parts = line.split()
            if len(parts) >= 6:
                rows.append(parts)
    return rows

fam = read_fam(bfile) if bfile else []
bim = read_bim(bfile) if bfile else []
n = len(fam)

if "-gk" in args:
    # 全队列 K：严格对称 PSD（对角 1 + 小值）
    k = [[1.0 if i == j else 0.05 for j in range(n)] for i in range(n)]
    suffix = "cXX" if "1" in args[args.index("-gk") + 1:] or args[args.index("-gk") + 1] == "1" else "sXX"
    with open(os.path.join(outdir, outname + f".{suffix}.txt"), "w") as h:
        for row in k:
            h.write("\t".join(f"{v:.6f}" for v in row) + "\n")
    with open(os.path.join(outdir, outname + ".log.txt"), "w") as h:
        h.write("mock gemma -gk done\n")
    sys.exit(0)

if "-lmm" in args:
    # assoc 输出（wald：.assoc.txt）
    lmm_code = "1"
    idx = args.index("-lmm")
    if idx + 1 < len(args):
        lmm_code = args[idx + 1]
    suffix_map = {"1": "assoc", "2": "lrt", "3": "score", "4": "assoc"}
    suffix = suffix_map.get(lmm_code, "assoc")
    with open(os.path.join(outdir, outname + f".{suffix}.txt"), "w") as h:
        h.write("chr\trs\tn_miss\tn_obs\tallele1\tallele2\taf\tbeta\tse"
                "\tlogl_H1\tl_remle\tp_wald\n")
        for r in bim:
            chrom, vid, _cm, pos, a1, a2 = r[0], r[1], r[2], r[3], r[4], r[5]
            af = rng.uniform(0.05, 0.5)
            beta = rng.gauss(0, 0.2)
            se = abs(rng.gauss(0.05, 0.01))
            p = rng.uniform(0.05, 1)
            h.write(f"{chrom}\t{vid}\t0\t{n}\t{a1}\t{a2}\t{af:.4f}"
                    f"\t{beta:.4f}\t{se:.4f}\t-100\t-100\t{p:.6g}\n")
    with open(os.path.join(outdir, outname + ".log.txt"), "w") as h:
        h.write("mock gemma -lmm done\n")
    sys.exit(0)

print("mock gemma: unsupported mode", file=sys.stderr)
sys.exit(1)
'''

BCFTOOLS_MOCK = r'''#!/usr/bin/env python3
# mock bcftools（GWASFlow L3 行为测试）
import gzip
import os
import sys
from pathlib import Path

args = sys.argv[1:]

def io_paths():
    src = dst = None
    for i, a in enumerate(args):
        if a in ("-o", "--output") and i + 1 < len(args):
            dst = args[i + 1]
    # 输入为最后一个非选项参数
    for a in reversed(args):
        if not a.startswith("-") and (Path(a).exists()) and a != dst:
            src = a
            break
    return src, dst

src, dst = io_paths()

sub = args[0] if args else ""

if sub == "view" and "-h" in args:
    # header 检查（DP/GQ 字段探测）
    has = False
    if src and (str(src).endswith(".gz") or str(src).endswith(".bcf")):
        try:
            with gzip.open(src, "rt", errors="replace") as h:
                for line in h:
                    if line.startswith("##FORMAT") and ("ID=DP" in line or "ID=GQ" in line):
                        has = True
                    if not line.startswith("#"):
                        break
        except OSError:
            pass
    if has:
        print("##FORMAT=<ID=DP,Number=1,Type=Integer>")
    sys.exit(0)

if sub == "query":
    fmt = ""
    for i, a in enumerate(args):
        if a == "-f" and i + 1 < len(args):
            fmt = args[i + 1]
    rows = []
    if src and Path(src).exists():
        opener = gzip.open if str(src).endswith(".gz") else open
        with opener(src, "rt", errors="replace") as h:
            for line in h:
                if line.startswith("#"):
                    continue
                parts = line.split("\t")
                chrom, pos, vid, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
                alts = alt.strip().split(",")
                if "%ID" in fmt:
                    rows.append(f"{vid}\t{chrom}\t{pos}\t{ref}\t{alts[0]}")
                elif "%CHROM" in fmt and "%REF" in fmt:
                    rows.append(f"{chrom}\t{ref}\t{alt.strip()}")
                elif "%CHROM" in fmt:
                    rows.append(chrom)
    sys.stdout.write("\n".join(rows) + ("\n" if rows else ""))
    sys.exit(0)

if sub in ("view", "norm", "annotate", "filter", "index"):
    import shutil as _shutil
    if dst:
        if src and Path(src).exists():
            rename_map_arg = None
            for i, a in enumerate(args):
                if a == "--rename-chrs" and i + 1 < len(args):
                    rename_map_arg = args[i + 1]
            set_id_fmt = None
            for i, a in enumerate(args):
                if a == "--set-id" and i + 1 < len(args):
                    set_id_fmt = args[i + 1]
            if set_id_fmt and src and Path(src).exists():
                # annotate --set-id：重写 ID 列（%CHROM:%POS:%REF:%FIRST_ALT）
                opener = gzip.open if str(src).endswith(".gz") else open
                out_lines = []
                with opener(src, "rt", errors="replace") as h:
                    for line in h:
                        if line.startswith("##"):
                            out_lines.append(line)
                        elif line.startswith("#CHROM"):
                            out_lines.append(line)
                        else:
                            parts = line.rstrip("\n").split("\t")
                            if len(parts) >= 5:
                                chrom, pos, _vid, ref, alt = parts[:5]
                                parts[2] = f"{chrom}:{pos}:{ref}:{alt.split(',')[0]}"
                                out_lines.append("\t".join(parts) + "\n")
                payload = "".join(out_lines).encode()
                if str(dst).endswith(".gz"):
                    with gzip.open(dst, "wb") as h:
                        h.write(payload)
                else:
                    Path(dst).write_bytes(payload)
            elif rename_map_arg and Path(rename_map_arg).exists():
                # 真实重命名（原染色体名 → 引擎数值名）
                mapping = {}
                for line in Path(rename_map_arg).read_text().splitlines():
                    parts = line.split("\t")
                    if len(parts) >= 2:
                        mapping[parts[0].strip()] = parts[1].strip()
                opener = gzip.open if str(src).endswith(".gz") else open
                out_lines = []
                with opener(src, "rt", errors="replace") as h:
                    for line in h:
                        if line.startswith("#"):
                            out_lines.append(line)
                        else:
                            parts = line.split("\t")
                            if parts and parts[0] in mapping:
                                parts[0] = mapping[parts[0]]
                            out_lines.append("\t".join(parts))
                payload = "".join(out_lines).encode()
                if str(dst).endswith(".gz"):
                    with gzip.open(dst, "wb") as h:
                        h.write(payload)
                else:
                    Path(dst).write_bytes(payload)
            else:
                # 其余子命令：字节级透传（输入已是 gzip；避免双层压缩）
                _shutil.copyfile(src, dst)
        else:
            Path(dst).write_bytes(b"")
        if str(dst).endswith(".gz") and sub != "index":
            Path(str(dst) + ".csi").touch()
        if sub == "index" and src:
            Path(str(src) + ".csi").touch()
    sys.exit(0)

print(f"mock bcftools: unsupported subcommand {sub}", file=sys.stderr)
sys.exit(1)
'''

RSCRIPT_MOCK = r'''#!/usr/bin/env python3
# mock Rscript（GWASFlow L3 行为测试）：生成最小合规图片/数据文件
import sys
from pathlib import Path

args = sys.argv[1:]
out_dir = None
for i, a in enumerate(args):
    if a == "--out-dir" and i + 1 < len(args):
        out_dir = args[i + 1]

PNG_MAGIC = bytes([0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A]) + b"\x00" * 8

if out_dir:
    d = Path(out_dir)
    d.mkdir(parents=True, exist_ok=True)
    for name in ["manhattan.svg", "manhattan.pdf", "manhattan.png",
                 "qq.svg", "qq.pdf", "qq.png"]:
        p = d / name
        if name.endswith(".png"):
            p.write_bytes(PNG_MAGIC)
        else:
            p.write_text(f"mock plot {name}\n")
    (d / "plot_data.tsv").write_text(
        "chrom\tpos\tvariant_id\tneg_log10_p_display\tp_underflow\tsignificant\n")
    # diagnostics 模式（无 gwas 输入时只生成 png）
    if not any(a == "--sumstats" for a in args):
        for name in ["pca.png", "kinship_heatmap.png", "phenotype_distribution.png",
                     "sample_inclusion.png"]:
            (d / name).write_bytes(PNG_MAGIC)
print("mock Rscript done")
sys.exit(0)
'''


def main(argv=None):
    parser = argparse.ArgumentParser(description="生成 GWASFlow L3 mock 工具")
    parser.add_argument("--bin-dir", required=True, help="mock bin 输出目录")
    args = parser.parse_args(argv)

    bin_dir = Path(args.bin_dir)
    bin_dir.mkdir(parents=True, exist_ok=True)
    for name, content in [
        ("plink2", PLINK2_MOCK),
        ("gemma", GEMMA_MOCK),
        ("bcftools", BCFTOOLS_MOCK),
        ("Rscript", RSCRIPT_MOCK),
    ]:
        path = bin_dir / name
        path.write_text(content)
        path.chmod(path.stat().st_mode | stat.S_IEXEC)
        print(f"[make_mock_tools] {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
