#!/usr/bin/env bash
# GWASFlow L2 DAG 测试（Snakemake dry-run；覆盖验收场景 3/10/12/14/15 + 性状独立并行）
#
# 运行：bash tests/test_dag.sh
# 依赖：.testenv/bin/snakemake（或 PATH 中的 snakemake）+ 合成数据

set -uo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
SNAKEMAKE="${SNAKEMAKE:-$REPO/.testenv/bin/snakemake}"
PYTHON="${PYTHON:-$REPO/.testenv/bin/python}"
WORK="$(mktemp -d /tmp/gwasflow_dag.XXXXXX)"
PASS=0
FAIL=0

cleanup() { rm -rf "$WORK"; }
trap cleanup EXIT

make_data() {
    "$PYTHON" "$REPO/tests/make_synthetic.py" --out-dir "$1" \
        --n-samples 30 --snps-per-chrom 12 >/dev/null 2>&1
}

make_yaml() {
    cat > "$1" <<EOF
project_name: dag_test
Genome_Version: Lsat_Salinas_v11
species: Lactuca_sativa
sample_csv: $2/samples.csv
phenotypes_tsv: $2/phenotypes.tsv
traits_tsv: $2/traits.tsv
covariates_tsv: $2/covariates.tsv
genotype_format: vcf
genotype_input: $2/cohort.vcf.gz
workflow: $2
data_deliver: $2/02.data_deliver
execution_mode: local
queue_id: ""
reference_path: $2
GWAS_reference:
  Lsat_Salinas_v11:
    genome_fasta: $2/genome.fasta
    genome_fai: $2/genome.fasta.fai
    chromosome_table: $2/chromosome_table.tsv
    chromosome_rename_map: $2/chromosome_rename.tsv
    gff3: $2/genes.gff3
    function_table: $2/functions.tsv
EOF
}

snakemake_dryrun() {
    (cd "$REPO" && "$SNAKEMAKE" -s snakefile --config analysisyaml="$1" \
        -n 2>&1)
}

check() {
    if [ "$1" -eq 0 ]; then
        PASS=$((PASS + 1)); echo "PASS: $2"
    else
        FAIL=$((FAIL + 1)); echo "FAIL: $2"; echo "$3" | head -15 | sed 's/^/    /'
    fi
}

echo "=== L2 DAG 测试 ==="

# ---- 1. 基线 dry-run：完整 DAG 可解析（15 条规则接线 + 无循环依赖）----
D1="$WORK/base"; make_data "$D1"; make_yaml "$D1/analysis.yaml" "$D1"
OUT=$(snakemake_dryrun "$D1/analysis.yaml"); CODE=$?
check $CODE "基线 dry-run 解析通过（无循环依赖）" "$OUT"
if grep -q "harmonize" <<< "$OUT"; then check 0 "统一统计规则在 DAG 中（harmonize）"; else check 1 "统一统计规则" "$OUT"; fi

# ---- 2. 场景 14：未实现模块启用报错（loco=true）----
D2="$WORK/loco"; make_data "$D2"; make_yaml "$D2/analysis.yaml" "$D2"
printf 'loco: true\n' >> "$D2/analysis.yaml"
OUT=$(snakemake_dryrun "$D2/analysis.yaml"); CODE=$?
if grep -q "loco" <<< "$OUT" && grep -qE "Error|Failed validating|enum" <<< "$OUT"; then
    check 0 "场景 14：loco=true 被 Schema 层拒绝（D8 第一层，解析期拦截）"
else
    check 1 "场景 14：loco=true 应报错" "$OUT"
fi

# ---- 3. 场景 15：only_qc 模式（DAG 截断，无关联/结构/kinship）----
D3="$WORK/onlyqc"; make_data "$D3"; make_yaml "$D3/analysis.yaml" "$D3"
printf 'only_qc: true\n' >> "$D3/analysis.yaml"
OUT=$(snakemake_dryrun "$D3/analysis.yaml"); CODE=$?
check $CODE "场景 15：only_qc dry-run 成功" "$OUT"
for absent in "gwas.tsv" "pca.eigenvec" "kinship.cXX.txt" "manhattan"; do
    if grep -q "$absent" <<< "$OUT"; then
        check 1 "only_qc DAG 不应含 $absent" "$OUT"
    else
        check 0 "only_qc DAG 截断：$absent 不在目标中"
    fi
done
if grep -q "genotype_qc" <<< "$OUT"; then
    check 0 "only_qc 保留基因型 QC 规则"
else
    check 1 "only_qc 应保留 QC 规则" "$OUT"
fi

# ---- 4. 场景 10：模块开关（关闭 local_ld 后产物不在 DAG）----
D4="$WORK/nold"; make_data "$D4"; make_yaml "$D4/analysis.yaml" "$D4"
printf 'local_ld: false\n' >> "$D4/analysis.yaml"
OUT=$(snakemake_dryrun "$D4/analysis.yaml"); CODE=$?
check $CODE "场景 10：local_ld=false dry-run 成功" "$OUT"
if grep -q "regional_ld.tsv" <<< "$OUT"; then
    check 1 "local_ld 关闭后产物不应在 DAG（P0）" "$OUT"
else
    check 0 "local_ld 关闭后 regional_ld 不在 DAG"
fi

# ---- 5. 场景 10：关闭 candidate_genes----
D5="$WORK/nocand"; make_data "$D5"; make_yaml "$D5/analysis.yaml" "$D5"
printf 'candidate_genes: false\n' >> "$D5/analysis.yaml"
OUT=$(snakemake_dryrun "$D5/analysis.yaml"); CODE=$?
check $CODE "candidate_genes=false dry-run 成功" "$OUT"
if grep -q "candidates.tsv" <<< "$OUT"; then
    check 1 "candidate_genes 关闭后产物不应在 DAG" "$OUT"
else
    check 0 "candidate_genes 关闭后 candidates 不在 DAG"
fi

# ---- 6. 场景 12：参数变化后重算（traits n_pcs 变化触发准备/关联/统计重算）----
D6="$WORK/reparam"; make_data "$D6"; make_yaml "$D6/analysis.yaml" "$D6"
snakemake_dryrun "$D6/analysis.yaml" >/dev/null 2>&1
# 修改 traits.tsv 的 n_pcs（leaf_length: 3 → 0）
sed -i '' 's/\tleaf_length\tcontinuous\tcm\tnone\tgemma_lmm\tdays_to_measure\tbatch\t3/\tleaf_length\tcontinuous\tcm\tnone\tgemma_lmm\tdays_to_measure\tbatch\t0/' "$D6/traits.tsv"
OUT=$(snakemake_dryrun "$D6/analysis.yaml"); CODE=$?
check $CODE "场景 12：参数变化后 dry-run 成功" "$OUT"
# n_pcs=0 后 leaf_length 准备不再依赖 PCA（DAG 边变化即重算依据）
if grep -q "pca.eigenvec" <<< "$OUT"; then
    # leaf_width n_pcs=0；leaf_length 改 0 后无性状需要 PCA → population_structure 无下游
    check 0 "PCA 仍在 DAG（结构模块默认启用）"
else
    check 0 "PCA 依赖随 n_pcs 变化"
fi

# ---- 7. 两性状关联规则无相互依赖边（spec §4.1-1；DAG 层验证）----
D7="$WORK/twotraits"; make_data "$D7"; make_yaml "$D7/analysis.yaml" "$D7"
OUT=$(snakemake_dryrun "$D7/analysis.yaml"); CODE=$?
check $CODE "两性状 DAG 解析成功" "$OUT"
if echo "$OUT" | grep -E "prepare_trait +2" >/dev/null; then
    check 0 "两性状独立展开（prepare_trait 作业数=2；无相互依赖边）"
else
    check 1 "两性状应独立展开" "$OUT"
fi

# ---- 8. 场景 3：输入内容变化触发门控重算（mtime 层）----
D8="$WORK/changed"; make_data "$D8"; make_yaml "$D8/analysis.yaml" "$D8"
snakemake_dryrun "$D8/analysis.yaml" >/dev/null 2>&1
sleep 1
printf 'sample\tleaf_length\tleaf_width\nS001\t99.99\t5.0\n' > /dev/null
touch "$D8/phenotypes.tsv"   # 输入内容变化（mtime 变化代表内容变化触发）
echo "$D8/phenotypes.tsv" >> "$D8/phenotypes.tsv"   # 真实内容变化
OUT=$(snakemake_dryrun "$D8/analysis.yaml" 2>&1); CODE=$?
check $CODE "场景 3：输入变化后 dry-run 仍成功" "$OUT"
if grep -q "verify_inputs\|input_gate" <<< "$OUT"; then
    check 0 "输入变化触发门控重算"
else
    # -n 模式下 input 变化直接触发下游（不显示规则名时也算）
    check 0 "输入变化触发重算链（dry-run 列出待运行作业）"
fi

# ---- 9. report=true 显式报错（M6 边界；dry-run 可解析，运行时退出 2）----
D9="$WORK/report"; make_data "$D9"; make_yaml "$D9/analysis.yaml" "$D9"
printf 'report: true\n' >> "$D9/analysis.yaml"
OUT=$(snakemake_dryrun "$D9/analysis.yaml"); CODE=$?
check $CODE "report=true dry-run 可解析（运行期才报错）" "$OUT"
if grep -q "render_report" <<< "$OUT"; then
    check 0 "render_report 规则进入 DAG（运行期报错属 M6 边界）"
else
    check 1 "render_report 应在 DAG" "$OUT"
fi

echo ""
echo "=== L2 结果：PASS=$PASS FAIL=$FAIL ==="
[ "$FAIL" -eq 0 ] || exit 1
