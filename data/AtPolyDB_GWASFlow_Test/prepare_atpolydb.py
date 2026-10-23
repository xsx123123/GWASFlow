#!/usr/bin/env python3
"""Convert a verified AtPolyDB TAIR9 snapshot to matched PLINK TPED test inputs."""
import argparse
import collections
import csv
import hashlib
import io
import json
import math
import sys
import tarfile
from pathlib import Path

MEMBER = "call_method_75/call_method_75_TAIR9.csv"
MISSING = {"", "NA", "NAN", "."}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_phenotypes(path, column):
    raw = Path(path).read_bytes()
    try:
        text = raw.decode("utf-8-sig")
        encoding = "utf-8-sig"
    except UnicodeDecodeError:
        text = raw.decode("latin-1")
        encoding = "latin-1"
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    if not {"ecotype id", column} <= set(reader.fieldnames or []):
        raise ValueError(f"Required columns: ecotype id, {column}; observed: {reader.fieldnames}")
    if len(reader.fieldnames) != len(set(reader.fieldnames)):
        raise ValueError("Duplicate phenotype column names")
    result = {}
    for row in reader:
        sid = row["ecotype id"].strip()
        if not sid or sid in result or None in row:
            raise ValueError(f"Malformed phenotype row or duplicate ecotype ID: {sid}")
        raw_value = row[column].strip()
        value = None if raw_value.upper() in MISSING else float(raw_value)
        if value is not None and not math.isfinite(value):
            raise ValueError(f"Nonfinite phenotype: {sid}")
        result[sid] = value
    return result, encoding


def write_tsv(path, header, rows):
    with Path(path).open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def prepare(args):
    if args.snp_step < 1:
        raise ValueError("--snp-step must be at least 1")
    if args.min_samples < 3:
        raise ValueError("--min-samples must be at least 3")
    out = Path(args.out_dir)
    if out.exists() and any(out.iterdir()):
        raise ValueError("--out-dir must be new or empty; use a separate directory for each run")
    phenotype, encoding = read_phenotypes(args.phenotypes, args.trait)
    trait_name = args.trait_name
    if not trait_name or any(c.isspace() for c in trait_name):
        raise ValueError("--trait-name must be nonempty without whitespace")
    seen_positions = set()
    counts = collections.Counter()
    chr_input = collections.Counter()
    chr_output = collections.Counter()
    tokens = collections.Counter()
    with tarfile.open(args.archive, mode="r|gz") as archive:
        for member in archive:
            if member.name != MEMBER:
                continue
            reader = csv.reader(line.decode("ascii") for line in archive.extractfile(member))
            array_row = next(reader)
            ecotype_row = next(reader)
            if array_row[:2] != ["-", "-"] or ecotype_row[:2] != ["Chromosome", "Positions"]:
                raise ValueError("Unexpected genotype header: expected array IDs followed by ecotype IDs")
            if len(array_row) != len(ecotype_row):
                raise ValueError("Genotype header widths differ")
            ids = ecotype_row[2:]
            if not ids or any(not sid or any(c.isspace() for c in sid) for sid in ids):
                raise ValueError("Invalid genotype ecotype IDs")
            duplicates = [sid for sid, n in collections.Counter(ids).items() if n > 1]
            if duplicates:
                raise ValueError(f"Repeated genotype ecotype IDs; resolve technical replicates explicitly: {duplicates}")
            matched = [sid for sid in ids if sid in phenotype]
            selected = [sid for sid in matched if phenotype[sid] is not None]
            if len(selected) < args.min_samples:
                raise ValueError(f"Only {len(selected)} matched nonmissing samples; minimum {args.min_samples}")
            if len({phenotype[sid] for sid in selected}) < 2:
                raise ValueError("Selected phenotype has zero variance")
            selected_set = set(selected)
            indexes = [i + 2 for i, sid in enumerate(ids) if sid in selected_set]
            out.mkdir(parents=True, exist_ok=True)
            prefix = out / "atpolydb"
            with prefix.with_suffix(".tfam").open("w") as fam:
                for sid in selected:
                    fam.write(f"0\t{sid}\t0\t0\t0\t-9\n")
            write_tsv(out / "samples.tsv", ["sample", "fid", "iid", "include"],
                      [(sid, "0", sid, "true") for sid in selected])
            write_tsv(out / "phenotypes.tsv", ["sample", trait_name],
                      [(sid, format(phenotype[sid], ".17g")) for sid in selected])
            write_tsv(out / "phenotype.plink.tsv", ["FID", "IID", trait_name],
                      [("0", sid, format(phenotype[sid], ".17g")) for sid in selected])
            with prefix.with_suffix(".tped").open("w") as tped:
                for index, row in enumerate(reader):
                    counts["source_snps"] += 1
                    if len(row) != len(ecotype_row):
                        raise ValueError(f"Genotype row {index + 3} has the wrong number of fields")
                    chrom, pos = row[:2]
                    if chrom not in {"1", "2", "3", "4", "5"} or not pos.isdigit() or int(pos) < 1:
                        raise ValueError(f"Invalid nuclear coordinate: {chrom}:{pos}")
                    position = (chrom, int(pos))
                    if position in seen_positions:
                        raise ValueError(f"Duplicate coordinate: {chrom}:{pos}")
                    seen_positions.add(position)
                    chr_input[chrom] += 1
                    if index % args.snp_step:
                        counts["skipped_thinning"] += 1
                        continue
                    calls = [row[i].strip().upper() for i in indexes]
                    tokens.update(calls)
                    invalid = set(calls) - {"A", "C", "G", "T"}
                    if invalid:
                        raise ValueError(f"Unsupported genotype tokens at {chrom}:{pos}: {sorted(invalid)}; do not guess an encoding")
                    if len(set(calls)) == 1:
                        counts["skipped_monomorphic"] += 1
                        continue
                    if len(set(calls)) > 2:
                        counts["skipped_multiallelic"] += 1
                        continue
                    pairs = "\t".join(f"{allele}\t{allele}" for allele in calls)
                    tped.write(f"{chrom}\t{chrom}:{pos}\t0\t{pos}\t{pairs}\n")
                    counts["written_snps"] += 1
                    chr_output[chrom] += 1
            if not counts["written_snps"]:
                raise ValueError("No polymorphic biallelic SNPs remained")
            if set(chr_output) != {"1", "2", "3", "4", "5"}:
                raise ValueError("Fixture does not cover all five chromosomes; reduce --snp-step")
            report = {
                "status": "prepared_tped_not_yet_plink_or_gwas",
                "assembly": "TAIR9", "archive_member": MEMBER,
                "phenotype_source_column": args.trait, "phenotype_output_column": trait_name,
                "phenotype_encoding": encoding, "snp_step": args.snp_step,
                "genotype_samples": len(ids), "phenotype_samples": len(phenotype),
                "matched_samples": len(matched), "selected_nonmissing_samples": len(selected),
                "unmatched_phenotype_ids": sorted(set(phenotype) - set(ids)),
                "matched_missing_phenotype_ids": [s for s in matched if phenotype[s] is None],
                "snps": dict(counts), "source_snps_per_chr": dict(chr_input),
                "written_snps_per_chr": dict(chr_output), "selected_genotype_tokens": dict(tokens),
                "id_policy": "Match second genotype header (ecotype ID), FID=0, IID=ecotype ID; retain genotype column order",
                "genotype_policy": "A/C/G/T calls represented as homozygous diploid pairs; unknown tokens fail",
                "allele_policy": "Observed nucleotide alleles; PLINK A1/A2 are not asserted to be REF/ALT",
                "sources": [
                    {"file": Path(args.archive).name, "sha256": sha256(args.archive)},
                    {"file": Path(args.phenotypes).name, "sha256": sha256(args.phenotypes)}],
            }
            report["outputs_sha256"] = {p.name: sha256(p) for p in sorted(out.iterdir()) if p.is_file()}
            (out / "preparation.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
            print(json.dumps({"samples": len(selected), "snps": counts["written_snps"],
                              "assembly": "TAIR9", "out_dir": str(out)}, ensure_ascii=False))
            return
    raise ValueError(f"Archive does not contain {MEMBER}; never substitute TAIR8 silently")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", required=True, help="Downloaded call_method_75.tar.gz")
    parser.add_argument("--phenotypes", required=True, help="Downloaded phenotype_published_raw.tsv")
    parser.add_argument("--trait", default="5_FT10", help="Exact original phenotype column")
    parser.add_argument("--trait-name", default="FT10", help="Output trait column; set explicitly when changing --trait")
    parser.add_argument("--snp-step", type=int, default=1, help="Keep every Nth source SNP, then exclude monomorphic/multiallelic sites")
    parser.add_argument("--min-samples", type=int, default=20)
    parser.add_argument("--out-dir", required=True, help="New or empty directory")
    args = parser.parse_args()
    try:
        prepare(args)
    except (ValueError, OSError, csv.Error, tarfile.TarError, StopIteration) as exc:
        print(f"ERROR: {exc}. An interrupted output directory is incomplete; use a fresh directory after fixing the issue.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
