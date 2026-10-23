# AtPolyDB 真实数据端到端测试指南（L4/L5）

版本：v0.1.0；日期：2026-10-23。
数据：拟南芥 AtPolyDB（TAIR9 坐标，call_method_75，SHA256SUMS 锁定快照）——**真实群体基因型与表型**，接替生菜占位数据完成 L4（外部软件对账）与 L5（端到端真实数据）验收。
基线提交：见 `git log`（v0.1.0 初始基线之上）；通用交接见 `docs/HANDOVER.md`（本文件是其 L4/L5 的数据具体化）。

---

## 1. 数据包事实（已核验，勿重复猜测）

| 项 | 值 | 来源 |
| --- | --- | --- |
| 物种/组装 | Arabidopsis thaliana，**TAIR9**（自交二倍体） | 数据包 README §1 |
| smoke 子集 | **194 样本（FID=0，IID=ecotype）、10,703 SNP、常染色体 1–5、每 20 个 SNP 取 1** | preparation.json |
| 性状 | **FT10**（10°C 开花时间，连续，单位天）；194 全部非缺失；**无协变量** | phenotypes.tsv |
| 现有格式 | TPED/TFAM（无显式等位基因列；纯合自交系编码） | atpolydb.tped |
| 等位基因语义 | BIM A1/A2 为**观测等位基因，无 REF/ALT 断言** | README §4（红线） |
| 快照校验 | 11 个文件 SHA256 全部 OK（2026-10-23 实测） | SHA256SUMS |
| 数据包合同测试 | 4/4 通过（ecotype 匹配/二倍体/TAIR8 拒绝/重复拒绝） | test_prepare_atpolydb.py |

**smoke 子集定位**：数据流/运行/绘图验证；**完整关联与局部 LD 验收建议用完整 SNP 集**（`--snp-step 1` 生成 `data/full`，见数据包 README §3）。抽稀子集不用于精细定位验收。

## 2. 与 GWASFlow 默认配置的关键差异（生菜→拟南芥）

| 配置项 | 生菜占位值 | AtPolyDB 真实值 | 说明 |
| --- | --- | --- | --- |
| `Genome_Version` | Lsat_Salinas_v11 | **TAIR9** | analysisyaml 覆盖 can_use_genome_version + GWAS_reference.TAIR9 |
| `parameter.qc.chr_count` | 9 | **5** | 拟南芥常染色体 1–5 |
| `genotype_format` | vcf | **bed** | 数据包产出 BED（TPED 经 PLINK 1.9 转换） |
| `candidate_genes` | true | **false** | 本包无 TAIR9 GFF3；诚实关闭。接入 TAIR9 注释后再启用，**不得混用 TAIR10 或生菜参考** |
| species | Lactuca_sativa | Arabidopsis_thaliana | — |

自交植物语义已由 GWASFlow 默认满足：HWE 关闭、高亲缘材料保留（拟南芥自交系预期大量高 K 值属正常，不是异常）。

## 3. 已完成的准备工作（本仓库内，可直接复用）

| 资产 | 位置 | 说明 |
| --- | --- | --- |
| 一键准备脚本 | `tests/atpolydb/prepare_inputs.sh` | samples.tsv→samples.csv + variant_map + PLINK 转换 + 占位协变量表 |
| variant_map 生成器 | `tests/atpolydb/make_variant_map_from_bim.py` | 支持 `--tped`（当前已实测：10,703 位点全登记）与 `--bim`（转换后交叉对账）双输入 |
| TAIR9 染色体表 | `tests/atpolydb/chromosome_table.tair9.tsv` | 绘图与区域展示（chr1–5 长度） |
| traits.tsv | `tests/atpolydb/traits.tsv` | FT10 九列契约（n_pcs=3；拟南芥 1001 群体结构著名，PCA 应出真实分层信号） |
| 配置模板 | `tests/atpolydb/analysisyaml.atpolydb_smoke.yaml` | 本机绝对路径，可直接 dry-run |
| 生成的 variant_map | `data/AtPolyDB_GWASFlow_Test/data/smoke/variant_map.tsv` | 观测等位基因登记（ref=字典序较小；集合与 BIM 对账一致） |

**接入时已修复的真实 bug**（真实数据暴露，已含代码修复+回归）：
1. BED 路线参考资源按需校验误要求 `chromosome_rename_map`（`rules/utils/reference_update.py` + `scripts/verify_inputs.py`）——BED 路线 BIM 已是数值染色体，已 discard。
2. schema `covariates_tsv` 不允许空值，无协变量项目（AtPolyDB 正是）无法配置——已允许 null。
3. samples 表契约差异：数据包 `samples.tsv` 为 TAB 分隔，GWASFlow 契约为逗号分隔 CSV——准备脚本自动转换。

## 4. 分步执行路径（在 HANDOVER §4 通用路径上落到 AtPolyDB 具体）

### Step 0：环境安装（当前机器 4 个外部软件均未装）

```bash
# Python/Snakemake（.testenv 已就绪；缺失则重建）
uv venv .testenv --python 3.12
uv pip install --python .testenv/bin/python pytest pyyaml "snakemake>=9.9.0"

# PLINK 1.9（数据包 TPED 转换用，命令名默认 plink）
conda install -c bioconda plink          # 或: conda create -n plink19 -c bioconda plink

# plink2 + gemma + bcftools（GWASFlow 生产链）
conda env create -f envs/plink2.yaml
conda env create -f envs/gemma.yaml
conda env create -f envs/bcftools.yaml

# R（绘图；smoke 全链含可视化时需要）
conda env create -f envs/r-plot.yaml

# 版本台账（安装后立即记录并回填 docs/VALIDATION.md §4）
plink --version | head -1; plink2 --version; gemma 2>&1 | head -1; bcftools --version | head -1; Rscript --version
```

### Step 1：输入准备与快照校验

```bash
cd data/AtPolyDB_GWASFlow_Test && shasum -a 256 -c SHA256SUMS && cd ../..
.testenv/bin/python -m unittest discover -s data/AtPolyDB_GWASFlow_Test -p "test_*.py" -v

bash tests/atpolydb/prepare_inputs.sh          # 核心准备（PLINK 装好后重跑补 BED 三件套）
```

就绪断言：samples.csv / phenotypes.tsv / variant_map.tsv / genotype.{bed,bim,fam} 全部 OK；
转换后核对 `variant_map.from_bim.tsv` 与 `variant_map.tsv` **位点集合一致**（交叉对账）。

### Step 2：dry-run（无软件依赖）

```bash
.testenv/bin/snakemake -s snakefile \
    --config analysisyaml=$PWD/tests/atpolydb/analysisyaml.atpolydb_smoke.yaml -n
```

断言：解析通过；生效模块开关清单中 `candidate_genes: disabled`（诚实关闭）；DAG 含
`normalize_genotype → genotype_qc → prune_pca → prepare_trait(FT10) → kinship → GEMMA → harmonize → loci/ld → plots → summary → deliver`。

### Step 3：only_qc 先行（门控+规范化+QC，不碰关联）

analysisyaml 追加 `only_qc: true` 后运行，核对：
- `01.gate/md5_check.tsv`：`baseline_only`（无 MD5 基线是预期；数据包完整性由 SHA256SUMS 独立保障）
- `02.genotype/normalize_audit.json`：`bed_route` 登记、`n_bim_variants=10703`、BIM/variant_map 对账通过
- `03.genotype_qc/qc_audit.json`：`hwe=off`、`--chr-set 5 no-x no-y no-xy no-mt`、家系字段不影响样本计数（FID=0 场景）
- 四阶段流失后样本数仍 ≥194×(1-mind) 量级（预期 near 194，自交系缺失率极低）

**FID=0 注意**：PLINK 语义中 FID=0 可能被视作缺失 FID；GWASFlow 已固定 `--nonfounders`
（纳入全部材料统计），若真实 plink2 输出样本数异常，优先检查该点（对应 HANDOVER §6 R10）。

### Step 4：全链运行（smoke）

```bash
.testenv/bin/snakemake -s snakefile \
    --config analysisyaml=$PWD/tests/atpolydb/analysisyaml.atpolydb_smoke.yaml \
    --cores 8 --use-conda --rerun-triggers mtime
```

重点核对（真实数据预期现象）：
- **PCA**：拟南芥 1001 群体有 9 大洲际群系结构，`12.plots/diagnostics/pca.png` 应见真实分层（合成数据没有的信号）
- **K 热图**：自交系块状高 K（对角+同群系块）属正常
- **FT10 关联**：FT10 已知主效位点区域（如 chr1 顶部 FLC 附近、chr4 FRI 区域在 TAIR9 坐标下）是否在 Manhattan 出现信号——**作生物学 sanity 参考，不作为唯一验收标准**（数据包 README §5 明示）
- `08.results/FT10/metadata.json`：M=有效 P 数（GEMMA -maf/-miss 过滤后，可能 <10703，unreported 表应可对账——验收场景 18）
- 无显著结果是合法科学结果（smoke 抽稀子集可能无 Bonferroni 显著位点）

### Step 5：L4 数值对账（验收场景 13 核心）

**口径 A（逐参数一致，数值必须逐位对上）**——复现 GWASFlow 实际命令（从 `run_metadata.json` 的 `command` 字段复制）：

```bash
# K（与 GWASFlow compute_kinship 同参数：-gk 1 -p 全1表型 -bfile cohort）
gemma -bfile <workflow>/03.genotype_qc/cohort -gk 1 \
    -p <workflow>/05.phenotype/FT10/kinship_all.pheno.txt -n 1 \
    -outdir /tmp/l4_ref -o ref_K

# LMM（与 run_gemma 同参数：-lmm 1(wald) -c 协变量(含截距+3PC) -maf 0.05 -miss 0.10 -hwe 0 -r2 0.9999）
gemma -bfile <workflow>/03.genotype_qc/cohort \
    -p <workflow>/05.phenotype/FT10/phenotype.gemma.txt \
    -c <workflow>/05.phenotype/FT10/covariates.gemma.txt \
    -k <workflow>/06.kinship/kinship.cXX.txt -lmm 1 \
    -maf 0.05 -miss 0.10 -hwe 0 -r2 0.9999 -n 1 \
    -outdir /tmp/l4_ref -o ref_assoc
```

断言：n_obs、beta、se、p_wald 与 `08.results/FT10/gwas.tsv` 数值一致（容差按版本声明，建议 1e-6 相对）；标记集合一致（cohort.bim 逐一比对）。

**口径 B（数据包 README 基线，仅 sanity 对比）**：`-lmm 4` 且**不加 PC 协变量**的基线命令（数据包 README §5）——与含 PC 的正式模型**不应要求 P 值完全一致**（README 自身声明）；只对比显著峰位置与 P 数量级。**防止误判**：口径 B 不一致不是 bug。

**口径 C（QC 对账）**：plink2 独立命令复现四阶段过滤（`tests/commands/README.md` 模板，chr_count=5），比对每阶段剔除位点/样本集合与 `removed_*.tsv` 一致。

### Step 6：完整集复跑（L5 加严）

```bash
# 生成完整 SNP 集（需下载原始数据，见数据包 README §2；约 92MB 下载 + 磁盘空间）
bash data/AtPolyDB_GWASFlow_Test/download_atpolydb.sh data/AtPolyDB_GWASFlow_Test/data/raw
.testenv/bin/python data/AtPolyDB_GWASFlow_Test/prepare_atpolydb.py \
    --archive data/AtPolyDB_GWASFlow_Test/data/raw/call_method_75.tar.gz \
    --phenotypes data/AtPolyDB_GWASFlow_Test/data/raw/phenotype_published_raw.tsv \
    --trait 5_FT10 --trait-name FT10 --snp-step 1 \
    --out-dir data/AtPolyDB_GWASFlow_Test/data/full
(cd data/AtPolyDB_GWASFlow_Test && bash make_plink.sh data/full)
bash tests/atpolydb/prepare_inputs.sh --full

# 复制 analysisyaml 改 genotype_bfile/输入路径指向 data/full 后全链重跑
```

完整集（~214k SNP）验证：LD 剪枝窗口语义、loci clumping 真实行为、大规模 BIM 的 variant_map 生成耗时、内存排序规模路线边界。

### Step 7：回归 + 记录

```bash
.testenv/bin/python -m pytest tests/test_unit.py -q   # L1
bash tests/test_dag.sh                                 # L2
bash tests/test_behavior.sh                            # L3（mock 层不受真实软件影响）
.testenv/bin/python -m unittest discover -s data/AtPolyDB_GWASFlow_Test -p "test_*.py" -v  # 数据包合同测试
```

更新 `docs/VALIDATION.md`：L4/L5 结果（分口径 A/B/C）、版本台账、未验证项收缩；新发现的问题与修复。

## 5. AtPolyDB 专用高风险区（补充 HANDOVER §6 R1–R10）

| # | 位置 | 风险 | 指引 |
| --- | --- | --- | --- |
| A1 | PLINK 1.9 tped 转换 | `--keep-allele-order` 下 A1/A2 分配规则未离线确认 | 转换后跑 `make_variant_map_from_bim.py --bim` 交叉对账：位点**集合**一致即通过（harmonize 对账是集合匹配，方向由 reported 语义保护） |
| A2 | FID=0 语义 | plink2 可能将 FID=0 视为缺失 FID，影响 founder 判定/警告 | GWASFlow 已固定 `--nonfounders`；若样本计数异常查 HANDOVER R10 |
| A3 | 自交系杂合缺失 | call_method_75 个别杂合 call 存在（tped 双列不同字母） | 正常现象；`prepare_trait_inputs` 与 QC 均按标准语义处理，无需干预 |
| A4 | GEMMA K 对角元 | 自交系近交系数高，centered K 可能有负特征值超出数值容差 | `validate_kinship` 只容忍数值误差级别；若 PSD 检查失败属**真实数据特性**而非流程 bug，记录并交方法审阅（禁止静默谱裁剪） |
| A5 | 基线口径混淆 | 数据包 README 基线是 `-lmm 4` 无 PC | 见 §4 Step 5 口径 B：不作逐位对账，防误判为 bug |
| A6 | FT10 生物学先验 | 已知主效位点（FLC/FRI 区域）可能出现在完整集 | 作 sanity 参考不作验收标准（数据包 README §5） |
| A7 | TAIR9 注释获取 | candidate_genes 当前诚实关闭 | 启用须配 TAIR9 GFF3（如 TAIR9_gene_models）；**不得混用 TAIR10 或生菜参考**（数据包 README §5 红线） |
| A8 | run 工作区位置 | workflow 若设在数据目录会污染入库快照 | 模板已用 `run/01.workflow`（.gitignore 已排除）；勿改为 `data/smoke/` 本身 |

## 6. 验收完成标准（AtPolyDB 口径的 L4/L5 定义）

在 HANDOVER §7 通用标准之上，AtPolyDB 特有：
1. 口径 A 对账逐位通过（N/beta/SE/P/标记集，容差记录进 VALIDATION.md）
2. 口径 C 四阶段 QC 剔除集合逐一对上
3. PCA 呈现真实群体结构（与拟南芥已知群系一致性作定性 sanity）
4. unreported_variants 与 GEMMA 过滤日志可对账（场景 18）
5. 完整集（Step 6）至少完成一次全链 + 口径 A 对账（smoke 仅作数据流验证）
6. 数据包 SHA256SUMS 快照校验记录在案；版本台账完整

## 7. 快速命令备忘

```bash
export ATY=$PWD/tests/atpolydb/analysisyaml.atpolydb_smoke.yaml

# dry-run / only_qc / 全链
.testenv/bin/snakemake -s snakefile --config analysisyaml=$ATY -n
.testenv/bin/snakemake -s snakefile --config analysisyaml=$ATY --cores 8 --use-conda

# 单目标调试
.testenv/bin/snakemake -s snakefile --config analysisyaml=$ATY 08.results/FT10/gwas.tsv

# 失败恢复（不清空既有产物）
.testenv/bin/snakemake -s snakefile --config analysisyaml=$ATY --cores 8 --rerun-incomplete

# 数据准备 / 快照校验 / 数据包合同测试
bash tests/atpolydb/prepare_inputs.sh
(cd data/AtPolyDB_GWASFlow_Test && shasum -a 256 -c SHA256SUMS)
.testenv/bin/python -m unittest discover -s data/AtPolyDB_GWASFlow_Test -p "test_*.py" -v
```
