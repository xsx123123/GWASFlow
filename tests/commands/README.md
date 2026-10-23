# L4 外部软件对账命令（M3 具备时执行）

本目录存放与生产流程**同参数**的独立对照命令（参数手册 §3/§4 契约），用于
GEMMA 数值对账（验收场景 13）与软件过滤对账（场景 18）。

## 执行条件

```bash
command -v plink2 && command -v gemma && command -v bcftools   # 三者齐备
plink2 --version; gemma 2>&1 | head -1; bcftools --version     # 固定版本记录
```

## 独立对照命令模板（与 run_gemma.py 命令数组逐参数一致）

### QC 对照（01_genotype_qc_reference.sh 语义）

```bash
plink2 --bfile 02.genotype/raw_engine \
    --geno 0.20 --chr-set 9 no-x no-y no-xy no-mt --nonfounders \
    --make-bed --out /tmp/ref_stage1
plink2 --bfile /tmp/ref_stage1 --mind 0.10 \
    --chr-set 9 no-x no-y no-xy no-mt --nonfounders \
    --make-bed --out /tmp/ref_stage2
plink2 --bfile /tmp/ref_stage2 --geno 0.10 --maf 0.05 --mac 5 \
    --chr-set 9 no-x no-y no-xy no-mt --nonfounders \
    --make-bed --out /tmp/ref_cohort
plink2 --bfile /tmp/ref_cohort --freq --out /tmp/ref_final
```

### GEMMA K 对照

```bash
gemma -bfile 03.genotype_qc/cohort -gk 1 \
    -p 05.phenotype/leaf_length/kinship_all.pheno.txt -n 1 \
    -outdir /tmp/ref_kinship -o ref_k
```

### GEMMA LMM 对照（与 run_metadata.json 的 command 字段逐参数比对）

```bash
gemma -bfile 03.genotype_qc/cohort \
    -p 05.phenotype/leaf_length/phenotype.gemma.txt \
    -c 05.phenotype/leaf_length/covariates.gemma.txt \
    -k 06.kinship/kinship.cXX.txt -lmm 1 \
    -maf 0.05 -miss 0.10 -hwe 0 -r2 0.9999 -n 1 \
    -outdir /tmp/ref_assoc -o ref_assoc
```

## 对账断言（核对 N/beta/SE/P/标记集）

- cohort.fam/bim 与独立命令输出逐一比对（N 与标记集合一致）
- assoc 的 n_obs/beta/se/p_wald 与 `08.results/{trait}/gwas.tsv` 数值一致
  （容差按软件版本声明确认）
- run_metadata.json 的 software_filters（-maf/-miss/-hwe/-r2）与命令行一致
- unreported_variants.tsv 与 GEMMA 日志过滤记录对账

**当前状态：未执行（M3：外部软件二进制与版本未定）——记录为未验证项，
不冒称已验证。**
