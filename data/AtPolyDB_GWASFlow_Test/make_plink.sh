#!/usr/bin/env bash
set -euo pipefail

input_dir="${1:-data/smoke}"
trait_name="${2:-FT10}"
plink_bin="${PLINK_BIN:-plink}"

"$plink_bin" \
  --tfile "$input_dir/atpolydb" \
  --chr-set 5 no-x no-y no-xy no-mt \
  --pheno "$input_dir/phenotype.plink.tsv" \
  --pheno-name "$trait_name" \
  --allow-no-sex \
  --keep-allele-order \
  --make-bed \
  --out "$input_dir/genotype"
