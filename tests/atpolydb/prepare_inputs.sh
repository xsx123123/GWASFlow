#!/usr/bin/env bash
# AtPolyDB（TAIR9）→ GWASFlow BED 路线输入一键准备
#
# 产出（data/smoke/ 下）：
#   samples.csv      GWASFlow 样本表契约（逗号分隔；源 samples.tsv 为 TAB）
#   variant_map.tsv  观测等位基因登记（tped 或 bim 推断；集合与 BIM 对账一致）
#   genotype.bed/bim/fam  需 PLINK 1.9（PLINK_BIN，默认 plink）；未安装则跳过并提示
#
# 用法：bash tests/atpolydb/prepare_inputs.sh [--full]
#   默认准备 data/smoke 快速子集；--full 使用 data/full（需先按数据包 README 生成）

set -euo pipefail

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
PYTHON="${PYTHON:-$REPO/.testenv/bin/python}"
PLINK_BIN="${PLINK_BIN:-plink}"
DATA_ROOT="$REPO/data/AtPolyDB_GWASFlow_Test"

SUBDIR="smoke"
if [ "${1:-}" = "--full" ]; then SUBDIR="full"; fi
SMOKE_DIR="$DATA_ROOT/data/$SUBDIR"

echo "=== AtPolyDB → GWASFlow 输入准备（$SUBDIR 子集） ==="

# 1. samples.tsv（TAB）→ samples.csv（GWASFlow 契约为逗号分隔 CSV）
if [ ! -f "$SMOKE_DIR/samples.csv" ]; then
    awk -F'\t' 'BEGIN{OFS=","} NF>0 {print $1,$2,$3,$4}' \
        "$SMOKE_DIR/samples.tsv" > "$SMOKE_DIR/samples.csv"
    echo "[prepare] samples.csv 生成: $(wc -l < "$SMOKE_DIR/samples.csv") 行"
else
    echo "[prepare] samples.csv 已存在，跳过"
fi

# 2. variant_map.tsv（tped 推断；plink 转换后可用 --bim 交叉对账）
if [ ! -f "$SMOKE_DIR/variant_map.tsv" ]; then
    "$PYTHON" "$REPO/tests/atpolydb/make_variant_map_from_bim.py" \
        --tped "$SMOKE_DIR/atpolydb.tped" \
        --chr-count 5 \
        --out "$SMOKE_DIR/variant_map.tsv"
else
    echo "[prepare] variant_map.tsv 已存在，跳过"
fi

# 3. PLINK 1.9 转换（TPED/TFAM → BED 三件套；表型填入 FAM 第 6 列，GWASFlow 不读该列）
if [ ! -f "$SMOKE_DIR/genotype.bim" ]; then
    if command -v "$PLINK_BIN" >/dev/null 2>&1; then
        (cd "$DATA_ROOT" && PLINK_BIN="$PLINK_BIN" bash make_plink.sh "data/$SUBDIR")
        echo "[prepare] genotype.bed/bim/fam 转换完成"
        # 转换后交叉对账（bim 模式与 tped 模式产物集合应一致）
        "$PYTHON" "$REPO/tests/atpolydb/make_variant_map_from_bim.py" \
            --bim "$SMOKE_DIR/genotype.bim" \
            --chr-count 5 \
            --out "$SMOKE_DIR/variant_map.from_bim.tsv"
        echo "[prepare] 交叉对账文件: variant_map.from_bim.tsv（集合应与 variant_map.tsv 一致）"
    else
        echo "[prepare] 跳过 PLINK 转换（未找到 ${PLINK_BIN}）。"
        echo "          安装 PLINK 1.9 后重跑本脚本，或手动执行："
        echo "          (cd $DATA_ROOT && bash make_plink.sh data/$SUBDIR)"
    fi
else
    echo "[prepare] genotype.bed/bim/fam 已存在，跳过"
fi

# 4. 无协变量占位表（GWASFlow covariates_tsv 配置为非空路径；性状未声明协变量时不被消费）
if [ ! -f "$SMOKE_DIR/covariates.tsv" ]; then
    printf 'sample\n' > "$SMOKE_DIR/covariates.tsv"
    echo "[prepare] covariates.tsv 占位表生成（仅表头；FT10 无协变量声明，不会被消费）"
fi

echo ""
echo "=== 就绪状态 ==="
for f in samples.csv phenotypes.tsv variant_map.tsv genotype.bim genotype.fam genotype.bed; do
    if [ -f "$SMOKE_DIR/$f" ]; then
        echo "  OK   ${f}"
    else
        echo "  MISS ${f}（见上方提示）"
    fi
done
echo ""
echo "配置模板: tests/atpolydb/analysisyaml.atpolydb_smoke.yaml"
echo "测试指南: docs/ATPOLYDB_TEST_GUIDE.md"
