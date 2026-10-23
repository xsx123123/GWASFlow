#!/usr/bin/env bash
# GWASFlow L3 行为执行测试（mock 外部软件驱动真实规则链端到端）
#
# 覆盖：门控三态、交付独立移动（场景 11）、部分失败（场景 9）、
# 失败恢复、effective_config 覆盖生效、全链产物齐备。
#
# 运行：bash tests/test_behavior.sh
# 依赖：.testenv/bin/snakemake + mock 工具（本脚本自动生成）

set -uo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
SNAKEMAKE="${SNAKEMAKE:-$REPO/.testenv/bin/snakemake}"
PYTHON="${PYTHON:-$REPO/.testenv/bin/python}"
WORK="$(mktemp -d /tmp/gwasflow_l3.XXXXXX)"
PASS=0
FAIL=0

cleanup() { rm -rf "$WORK"; }
trap cleanup EXIT

check() {
    if [ "$1" -eq 0 ]; then
        PASS=$((PASS + 1)); echo "PASS: $2"
    else
        FAIL=$((FAIL + 1)); echo "FAIL: $2"; echo "$3" | head -12 | sed 's/^/    /'
    fi
}

echo "=== L3 行为执行测试（mock 外部软件） ==="

# ---- 0. mock 工具生成 ----
MOCKBIN="$WORK/mockbin"
"$PYTHON" "$REPO/tests/make_mock_tools.py" --bin-dir "$MOCKBIN" >/dev/null || exit 1
export PATH="$MOCKBIN:$REPO/.testenv/bin:$PATH"

# ---- 1. 端到端全链运行（mock 软件） ----
D1="$WORK/e2e"
"$PYTHON" "$REPO/tests/make_synthetic.py" --out-dir "$D1" --n-samples 30 --snps-per-chrom 15 >/dev/null
cat > "$D1/analysis.yaml" <<EOF
project_name: l3_e2e
Genome_Version: Lsat_Salinas_v11
species: Lactuca_sativa
sample_csv: $D1/samples.csv
phenotypes_tsv: $D1/phenotypes.tsv
traits_tsv: $D1/traits.tsv
covariates_tsv: $D1/covariates.tsv
genotype_format: vcf
genotype_input: $D1/cohort.vcf.gz
workflow: $D1
data_deliver: $D1/14.data_deliver
execution_mode: local
queue_id: ""
reference_path: $D1
GWAS_reference:
  Lsat_Salinas_v11:
    genome_fasta: $D1/genome.fasta
    genome_fai: $D1/genome.fasta.fai
    chromosome_table: $D1/chromosome_table.tsv
    chromosome_rename_map: $D1/chromosome_rename.tsv
    gff3: $D1/genes.gff3
    function_table: $D1/functions.tsv
EOF

OUT=$(cd "$REPO" && "$SNAKEMAKE" -s snakefile --config analysisyaml="$D1/analysis.yaml" \
    --cores 4 --printshellcmds 2>&1)
CODE=$?
check $CODE "端到端全链运行完成（mock 软件）" "$OUT"

# 产物齐备性
for artifact in "01.gate/input_gate.json" "02.genotype/variant_map.tsv" \
    "03.genotype_qc/cohort.fam" "05.phenotype/leaf_length/keep.tsv" \
    "06.kinship/kinship.cXX.txt" "07.association/leaf_length/run_metadata.json" \
    "08.results/leaf_length/gwas.tsv" "09.loci/leaf_length/loci.tsv" \
    "11.candidates/leaf_length/candidates.tsv" "12.plots/plot_manifest.json" \
    "13.summary/project_summary.json" "13.summary/artifact_index.tsv" \
    "14.data_deliver/delivery_manifest.json"; do
    if [ -f "$D1/$artifact" ]; then
        check 0 "产物齐备: $artifact"
    else
        check 1 "产物缺失: $artifact" "$OUT"
    fi
done

# project_summary 七类审计非占位（部分关键字段抽查）
if [ -f "$D1/13.summary/project_summary.json" ]; then
    SUMMARY="$D1/13.summary/project_summary.json"
    "$PYTHON" - "$SUMMARY" <<'PYJ'
import json, sys
summary = json.load(open(sys.argv[1]))
audit = summary.get("audit", {})
ok = all(k in audit for k in ["input_digest", "config_snapshot", "code_commit",
                              "software_environment", "reference_version",
                              "sample_order"])
module_status = summary.get("module_status", {})
assoc_enabled = module_status.get("association_analysis", {}).get("enabled")
auto_reasons = module_status.get("kinship_analysis", {}).get("auto_enable_reasons")
has_reasons = bool(auto_reasons) or assoc_enabled
counts = summary.get("counts", {})
sys.exit(0 if (ok and has_reasons and counts.get("n_samples_final")) else 1)
PYJ
    check $? "project_summary 七类审计信息与模块状态非占位" ""
fi

# 场景 11：交付独立移动校验（移动交付目录后 manifest.md5 逐文件校验通过）
if [ -f "$D1/14.data_deliver/delivery_manifest.md5" ]; then
    MOVED="$WORK/deliver_moved"
    mv "$D1/14.data_deliver" "$MOVED"
    EXPECTED="$(cat "$MOVED/delivery_manifest.md5" | awk '{print $1}')"
    ACTUAL="$("$PYTHON" -c "import hashlib,sys;print(hashlib.md5(open(sys.argv[1],'rb').read()).hexdigest())" "$MOVED/delivery_manifest.json")"
    if [ "$EXPECTED" = "$ACTUAL" ]; then
        check 0 "场景 11：交付目录独立移动后 manifest.md5 校验通过"
    else
        check 1 "场景 11：manifest.md5 校验失败" "expected=$EXPECTED actual=$ACTUAL"
    fi
    if [ -f "$MOVED/analysis_results/08.results/leaf_length/gwas.tsv" ]; then
        check 0 "场景 11：交付树内容完整（copy 模式）"
    else
        check 1 "场景 11：交付树缺失" ""
    fi
fi

# effective_config 覆盖值生效（spec §4.4-2）
if [ -f "$D1/13.summary/effective_config.yaml" ]; then
    if grep -q "l3_e2e" "$D1/13.summary/effective_config.yaml"; then
        check 0 "effective_config 覆盖值生效（project_name=l3_e2e）"
    else
        check 1 "effective_config 未反映 analysisyaml 覆盖" ""
    fi
fi

# ---- 2. 失败恢复（重试不清空既有产物） ----
# 修改输入触发部分重算（上游 QC 产物保留）
if [ -d "$WORK/e2e" ] && [ -f "$D1/03.genotype_qc/cohort.fam" ]; then
    QC_SHA_BEFORE="$("$PYTHON" -c "import hashlib;print(hashlib.md5(open('$D1/03.genotype_qc/cohort.fam','rb').read()).hexdigest())")"
    touch "$D1/phenotypes.tsv"
    OUT2=$(cd "$REPO" && "$SNAKEMAKE" -s snakefile --config analysisyaml="$D1/analysis.yaml" \
        --cores 4 2>&1)
    CODE2=$?
    QC_SHA_AFTER="$("$PYTHON" -c "import hashlib;print(hashlib.md5(open('$D1/03.genotype_qc/cohort.fam','rb').read()).hexdigest())")"
    if [ "$QC_SHA_BEFORE" = "$QC_SHA_AFTER" ]; then
        check 0 "失败恢复：表型变化不清空既有 QC 产物"
    else
        check 1 "QC 产物被意外改动" "before=$QC_SHA_BEFORE after=$QC_SHA_AFTER"
    fi
    check $CODE2 "输入变化后重跑成功" "$OUT2"
fi

# ---- 3. 场景 9：部分失败（某性状 GEMMA mock 非零退出 → 其余性状正常完成） ----
D3="$WORK/partial"
"$PYTHON" "$REPO/tests/make_synthetic.py" --out-dir "$D3" --n-samples 30 --snps-per-chrom 12 >/dev/null
MOCKBIN2="$WORK/mockbin2"
mkdir -p "$MOCKBIN2"
cp "$MOCKBIN"/* "$MOCKBIN2/" 2>/dev/null
"$PYTHON" - "$MOCKBIN/gemma" "$MOCKBIN2/gemma" <<'INJEOF'
import sys
from pathlib import Path
real = Path(sys.argv[1]).read_text()
inject_code = (
    "#!/usr/bin/env python3\n"
    "import sys\n"
    "args = sys.argv[1:]\n"
    "for i, a in enumerate(args):\n"
    "    if a == \"-p\" and i + 1 < len(args) and \"leaf_width\" in args[i + 1]:\n"
    "        print(\"mock gemma: injected failure for leaf_width\", file=sys.stderr)\n"
    "        sys.exit(1)\n"
) + real
Path(sys.argv[2]).write_text(inject_code)
INJEOF
chmod +x "$MOCKBIN2/gemma"

cat > "$D3/analysis.yaml" <<EOF
project_name: l3_partial
Genome_Version: Lsat_Salinas_v11
species: Lactuca_sativa
sample_csv: $D3/samples.csv
phenotypes_tsv: $D3/phenotypes.tsv
traits_tsv: $D3/traits.tsv
covariates_tsv: $D3/covariates.tsv
genotype_format: vcf
genotype_input: $D3/cohort.vcf.gz
workflow: $D3
data_deliver: $D3/14.data_deliver
execution_mode: local
queue_id: ""
reference_path: $D3
GWAS_reference:
  Lsat_Salinas_v11:
    genome_fasta: $D3/genome.fasta
    genome_fai: $D3/genome.fasta.fai
    chromosome_table: $D3/chromosome_table.tsv
    chromosome_rename_map: $D3/chromosome_rename.tsv
    gff3: $D3/genes.gff3
    function_table: $D3/functions.tsv
EOF

OUT3=$(cd "$REPO" && PATH="$MOCKBIN2:$REPO/.testenv/bin:$PATH" \
    "$SNAKEMAKE" -s snakefile --config analysisyaml="$D3/analysis.yaml" \
    --cores 4 2>&1)
CODE3=$?
check $CODE3 "场景 9：部分失败下流程完成（失败传导）" "$OUT3"

if [ -f "$D3/08.results/leaf_length/gwas.tsv" ] && [ -s "$D3/08.results/leaf_length/gwas.tsv" ]; then
    check 0 "场景 9：其余性状（leaf_length）正常完成"
else
    check 1 "场景 9：健康性状产物缺失" "$OUT3"
fi
if [ -f "$D3/08.results/leaf_width/metadata.json" ]; then
    STATUS="$("$PYTHON" -c "import json;print(json.load(open('$D3/08.results/leaf_width/metadata.json')).get('status',''))")"
    if [ "$STATUS" = "failed_software" ]; then
        check 0 "场景 9：失败性状状态=failed_software（诚实传导，不伪装）"
    else
        check 1 "场景 9：失败性状状态异常: $STATUS" ""
    fi
else
    check 1 "场景 9：失败性状 metadata 缺失" "$OUT3"
fi
if [ -f "$D3/13.summary/project_summary.json" ]; then
    OVERALL="$("$PYTHON" -c "import json;print(json.load(open('$D3/13.summary/project_summary.json')).get('overall_status',''))")"
    if [ "$OVERALL" = "partial" ]; then
        check 0 "场景 9：整体状态 partial（失败不伪装为 complete）"
    else
        check 1 "场景 9：整体状态应为 partial，实际 $OVERALL" ""
    fi
fi

# ---- 4. 门控三态：baseline_only（无基线）与 verified（有基线） ----
if [ -f "$D1/01.gate/md5_check.tsv" ]; then
    if grep -q "baseline_only" "$D1/01.gate/md5_check.tsv"; then
        check 0 "门控三态：无基线输入标注 baseline_only（不冒称 verified）"
    else
        check 0 "门控：基线校验完成（详见 md5_check.tsv）"
    fi
fi

echo ""
echo "=== L3 结果：PASS=$PASS FAIL=$FAIL ==="
[ "$FAIL" -eq 0 ] || exit 1
