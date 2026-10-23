# GWASFlow 交接文档：真实 GWAS 数据测试指南

版本：v0.1.0-handover；日期：2026-10-23。
交接目的：下一个执行 agent 使用**真实植物群体 GWAS 数据**完成 L4（外部软件对账）与 L5（端到端真实数据）验收。
基线提交：`afe633c`（工作区干净，78 个文件，全部 L1/L2/L3 测试通过）。

---

## 0. 必读顺序（执行 agent 先读这三份再动手）

1. `docs/VALIDATION.md` —— 已验证/未验证项的权威状态（L4/L5 全部未执行，不得冒称已验证）
2. `docs/SCRIPT_PARAMETERS.md` —— 全部脚本 CLI 契约与软件命令语义（参数逐字一致的红线）
3. `docs/GWASFlow_IMPLEMENTATION.md` —— 分析链路与方法约束（§9 质控约束、§10 GEMMA 主模型、§11 多重检验）

规格链（背景理解用，非执行必需）：`.codeartsdoer/specs/gwas_gemma_loop/{spec,design,tasks}.md`。

## 1. 当前基线状态

| 项 | 状态 |
| --- | --- |
| 代码 | v0.1.0 GEMMA 主分析闭环全链 15 条规则（提交 `afe633c`） |
| L1 单元（pytest） | 44/44 通过 |
| L2 DAG（snakemake dry-run） | 21/21 通过 |
| L3 mock 端到端 | 25/25 通过（**mock 软件**：`tests/make_mock_tools.py` 生成） |
| L4 真实软件对账 | **未执行**（本任务核心：验收场景 13 GEMMA 数值对账、18 软件过滤对账） |
| L5 真实数据端到端 | **未执行**（本任务核心） |
| 本机真实软件 | plink2 / gemma / bcftools / Rscript **均未安装**（`command -v` 检测） |

**红线**：真实数据测试的结论只能来自真实软件运行；mock 与合成数据（`tests/make_synthetic.py`）只能用于回归测试，不得写入 L4/L5 验收结论。

## 2. 运行环境准备（真实软件安装）

```bash
# 2.1 Python/Snakemake 测试环境（已存在于 .testenv/，若缺失则重建）
uv venv .testenv --python 3.12
uv pip install --python .testenv/bin/python pytest pyyaml "snakemake>=9.9.0"

# 2.2 外部软件 conda 环境（envs/ 已就位；版本占位——安装后立即固定并记录）
conda env create -f envs/gwasflow-genotype.yaml   # python3 + bcftools + plink2
conda env create -f envs/gemma.yaml               # python3 + gemma + openblas + numpy
conda env create -f envs/r-plot.yaml              # R + ggplot2 + jsonlite

# 2.3 版本固定（M3 消解的第一步：写入下方台账）
plink2 --version
gemma 2>&1 | head -1            # gemma 无 --version，无参数运行打印版本
bcftools --version | head -1
Rscript --version
```

**版本台账（agent 安装后填写并提交到 docs/VALIDATION.md §4 替换占位）**：

| 软件 | 安装版本 | 二进制路径 | 验收日期 |
| --- | --- | --- | --- |
| plink2 | 待填 | 待填 | 待填 |
| gemma | 待填 | 待填 | 待填 |
| bcftools | 待填 | 待填 | 待填 |
| R | 待填 | 待填 | 待填 |

## 3. 真实数据准备契约（用户提供方责任，逐项核对）

> 当前已提供的真实数据：`data/AtPolyDB_GWASFlow_Test/`（拟南芥 TAIR9）。该数据的全部准备
> 资产与接入步骤已就绪，**见 `docs/ATPOLYDB_TEST_GUIDE.md` §1–§3**；以下通用契约适用于
> 后续其他数据集接入。

### 3.1 基因型（二选一）

**VCF 路线（主路线）**：
- `cohort.vcf.gz`（+ `.csi` 或 `.tbi` 索引；**必须同组装下联合分型**，不得拼接单样本 VCF）
- 硬分型双等位 SNP（剂量输入会被 `--representation hardcall` 显式拒绝）
- 多等位/INDEL 会被排除并计数（`normalize_audit.json`）

**BED 路线**：
- 三件套 + **已核验的 variant_map.tsv**（`--input-variant-map`；BIM 的 A1/A2 不能证明 REF/ALT）

### 3.2 表格（全部 TSV/CSV，sample 唯一）

| 文件 | 必填列 | 说明 |
| --- | --- | --- |
| `samples.csv` | sample,fid,iid | VCF 路线 `--double-id` 要求 FID=IID=VCF 样本名；可选 include/group |
| `phenotypes.tsv` | sample + 每性状一列 | 缺失写 NA/NaN/.；**-9 是合法负值**（需视为缺失才加入 missing-values） |
| `traits.tsv` | 九列契约 | 见 `examples/tables/traits.tsv`；脚本参数由该表单点派生 |
| `covariates.tsv` | sample + 协变量列 | 数值/分类列必须显式区分声明于 traits.tsv |
| `expected_md5.tsv`（可选） | file/md5 | 提供方校验和基线；无基线时门控标注 baseline_only |

### 3.3 参考资源（按 Genome_Version 登记，绝对路径）

需准备并登记到 analysisyaml 的 `GWAS_reference.<version>`：

| 资源 | 用途 | 必需性 |
| --- | --- | --- |
| genome_fasta + fai | REF 校验（bcftools norm --check-ref e 严格拒绝不一致） | VCF 路线必需 |
| chromosome_rename.tsv | 原染色体名→引擎数值名（两列无表头；**必须覆盖输入 VCF 全部染色体**，未知染色体报错不压 0） | VCF 路线必需 |
| chromosome_table.tsv | 染色体顺序+长度（绘图与区域展示） | visualization 必需 |
| gff3 | 候选基因注释（gene feature + ID 属性） | candidate_genes 必需 |
| function_table | 基因功能描述（gene_id/description；一基因多行会汇总） | 可选 |

**M2 必须先解决**：生菜真实组装版本（当前占位 `Lsat_Salinas_v11`、9 对染色体）——真实数据到位后核对：
1. `Genome_Version` 与 reference.yaml `can_use_genome_version` 匹配
2. `parameter.qc.chr_count` = 真实染色体对数
3. rename map 覆盖 VCF 实际出现的全部染色体名（`bcftools query -f '%CHROM\n' cohort.vcf.gz | sort -u` 核对）

### 3.4 性状设计红线（不可由 agent 事后更改）

- 主检验预声明（`parameter.association.primary_test`，默认 wald）；**不按显著位点最多选择模型**
- `family_size` 预声明；部分性状失败不得事后缩小
- 有效材料 < min_samples（默认 20）阻断；该阈值不代表足够统计功效
- transform 仅 none/log/log1p（log 要求 y>0、log1p 要求 y>-1）

## 4. 分步执行路径（建议顺序，每步通过再进下一步）

### Step 0：干跑验证（无软件依赖）

```bash
cd /Users/zj/zj_code/GWASFlow
snakemake -s snakefile --config analysisyaml=/path/to/analysis.yaml -n
```

断言：解析通过；控制台打印生效模块开关清单（15 项）；无循环依赖错误。

### Step 1：only_qc 模式（先验门控+规范化+QC 链，不碰关联分析）

analysisyaml 追加 `only_qc: true`，运行后核对：
- `01.gate/md5_check.tsv` 三态正确（有基线=verified；无=baseline_only）
- `02.genotype/normalize_audit.json`：多等位/INDEL 排除计数合理、染色体映射全覆盖
- `03.genotype_qc/qc_audit.json`：`hwe=off`、四阶段流失明细合理（样本/标记数与预期同量级）
- `03.genotype_qc/cohort.fam` 的 N 合理（QC 后 ≥ min_samples_after_qc）

**此处极易因真实数据触发阻断**（见 §6 高风险区 R1/R2）。

### Step 2：结构+表型准备（仍是 only_qc=false 的前置）

```bash
snakemake -s snakefile --config analysisyaml=... 04.structure/pca.eigenvec
snakemake -s snakefile --config analysisyaml=... 05.phenotype/<trait>/metadata.json
```

核对：`pca.eigenvec` 表头含 #FID/IID；prepare 的 `sample_inclusion.tsv` 排除原因合理、`keep.tsv` 数量 = 有效 N。

### Step 3：kinship + 单性状 GEMMA + 统一统计

```bash
snakemake -s snakefile --config analysisyaml=... 06.kinship/kinship_validation.json
snakemake -s snakefile --config analysisyaml=... 08.results/<trait>/gwas.tsv
```

核对（对应验收场景 13/18 的前置）：
- `kinship_validation.json` 五项检查通过（维度/对称/身份/有限/PSD；PSD 只容忍数值误差级别）
- `run_metadata.json`：status=software_exited_zero、software_filters 四值（maf/miss/hwe=0/r2=0.9999）、五类绑定
- `08.results/<trait>/metadata.json`：M=有效 P 数、阈值=alpha/(family_size×M)、unreported 数可解释

### Step 4：L4 数值对账（验收场景 13 核心——必须执行）

按 `tests/commands/README.md` 的独立对照命令（与 run_metadata.json 的 `command` 字段逐参数一致）运行同版本 GEMMA，断言：
- N（n_obs）、beta、SE、P（p_wald）与 `gwas.tsv` 数值一致（容差按版本声明确认）
- 标记集合一致（cohort.bim 与独立命令输出逐一比对）
- GEMMA 自身过滤（-maf/-miss）造成的未报告位点与 `unreported_variants.tsv` 对账（场景 18）

### Step 5：全链 + 交付 + 多性状

```bash
snakemake -s snakefile --config analysisyaml=... --cores 8 --use-conda
```

断言：`13.summary/project_summary.json` 七类审计非占位、overall_status 如实（complete 或 partial）；移动 `14.data_deliver/` 后 `delivery_manifest.md5` 校验通过（场景 11 真实复验）。

### Step 6：回归 + 记录

```bash
.testenv/bin/python -m pytest tests/test_unit.py -q   # L1 不回归
bash tests/test_dag.sh                                 # L2 不回归
bash tests/test_behavior.sh                            # L3 不回归（mock 层不受真实软件影响）
```

更新 `docs/VALIDATION.md`：L4/L5 结果、版本台账、未验证项清单收缩；新增发现的问题与修复。

## 5. 配置模板（真实数据 analysisyaml 最小集）

```yaml
project_name: <真实项目名>
Genome_Version: <真实组装版本>        # M2：替换占位
species: <物种>
sample_csv: /abs/path/samples.csv     # 全部绝对路径
phenotypes_tsv: /abs/path/phenotypes.tsv
traits_tsv: /abs/path/traits.tsv
covariates_tsv: /abs/path/covariates.tsv
genotype_format: vcf                  # vcf | bcf | bed
genotype_input: /abs/path/cohort.vcf.gz
workflow: /abs/path/01.workflow       # 工作区（分析产物目录）
data_deliver: /abs/path/02.data_deliver
execution_mode: local
queue_id: ""
reference_path: /abs/path/reference_root
GWAS_reference:
  <真实组装版本>:
    genome_fasta: /abs/path/genome.fasta
    genome_fai: /abs/path/genome.fasta.fai
    chromosome_table: /abs/path/chromosome_table.tsv
    chromosome_rename_map: /abs/path/chromosome_rename.tsv
    gff3: /abs/path/genes.gff3
    function_table: /abs/path/functions.tsv   # 可选
parameter:
  qc:
    chr_count: <真实染色体对数>       # M2：替换占位 9
```

## 6. 高风险区：真实软件下预期需要修正的点（mock 无法覆盖）

**这是交接最重要的一节**。L3 mock 故意简化了真实软件行为，以下位置首次接真实软件时最可能暴露问题，按暴露概率排序：

| # | 位置 | 风险 | 排查/修正指引 |
| --- | --- | --- | --- |
| R1 | `scripts/normalize_genotype.py` bcftools 命令链 | `--apply-filters PASS`（pass_only 默认）会把 FILTER="." 的真实记录全部过滤掉；真实 VCF 常见 "." | 确认提供方 FILTER 语义；若数据以 "." 为主，analysisyaml/参数设 `filter_policy: pass_or_unfiltered`（语义：保留未评估记录，不称已通过质量过滤） |
| R2 | 同上 `--check-ref e` | 真实 REF 与 FASTA 不一致（组装错版/链向）会整链阻断（这是**设计行为**，不是 bug） | 错误日志含位点与双方碱基；先人工核对少量不一致位点确认组装版本（M2），禁止改用 `-w`/自动换碱基绕过 |
| R3 | `scripts/genotype_qc.py` `parse_plink_log_reasons` | mock 日志格式与真实 plink2 log 不同，原因行解析可能为空 | 原因表是补充信息（确定性集合对比才是主断言：removed_*.tsv 由 BIM/FAM 差集计算，不依赖日志）；若日志格式漂移只需放宽解析正则 |
| R4 | `scripts/run_gemma.py` GEMMA 输出命名 | 真实 GEMMA 输出文件名为 `<prefix>.assoc.txt` 与 `cXX/sXX` 后缀已按手册实现，但不同编译版可能差异 | 运行后若报"GEMMA 未产出 K 矩阵/assoc"，检查 `$outdir` 实际文件名并调整 `compute_kinship.py` 的 candidates 匹配 |
| R5 | `scripts/prepare_trait_inputs.py` PCA 匹配 | plink2 `--pca scols=fid` 真实输出表头/精度与 mock 不同 | 单元测试已覆盖 #FID/IID 解析；若真实输出含额外列需调整 `read_pcs` 的 pc_start 逻辑 |
| R6 | `scripts/define_loci.py` `.clumped` 解析 | `parse_clumped` 按简化列位解析，真实 plink2 .clumped 列更多（SP2 含 =value 分号列表） | 已按 SNP/SP2 索引解析；真实格式验证后修正列索引；解析失败会回退固定窗口并记录 fallback_reason（诚实降级不阻断） |
| R7 | `scripts/run_regional_ld.py` plink2 `--r2` 输出 | 真实 `.ld` 为.gz 或含注释行差异 | 若 `.ld` 未生成，检查 plink2 版本是否需要 `--r2 gz n-header` 等修饰 |
| R8 | `scripts/plot_gwas.R` | 真实 R 环境缺 ggplot2/jsonlite，或 dplyr 语义差异（代码仅用 ggplot2+jsonlite） | envs/r-plot.yaml 已声明依赖；`aes_string` 在新版 ggplot2 已弃用（警告不报错）；若报错改用 `.data[[pc_cols[1]]]` 语法 |
| R9 | `rules/12.visualization.smk` 染色体表 | 真实染色体名非 Lsat* 时 chromosome_table 顺序驱动绘图偏移 | 表第一列必须与 gwas.tsv 的 chrom 值**逐字一致**（含大小写） |
| R10 | plink2 `--chr-set` 语义 | 染色体数 >9 时需同步 chr_count；部分植物性染色体命名特殊 | 已固定 no-x no-y no-xy no-mt；真实性染色体数据不在本轮范围（人类/特殊 profile 属后续工作包） |

**修 bug 的边界**：修复必须在生产代码（scripts/ 或 rules/）内做且同步补测试；**禁止**为让流程"跑通"而：改统计红线（P 值/阈值/M 计算）、静默跳过身份对账、把 FAIL 状态写成 PASS、删 --check-ref 严格模式。

## 7. 验收通过标准（L4/L5 完成的定义）

1. 真实软件全链在真实数据上 exit 0（或如实的 partial——部分性状失败属合法结果，只要状态诚实）
2. 场景 13：GEMMA 独立命令对账 N/beta/SE/P/标记集一致（容差写入 VALIDATION.md）
3. 场景 18：run_metadata 软件四过滤生效值与命令行一致；unreported 与日志可对账
4. 无显著结果是**合法科学结果**（空表正常交付，不是失败）
5. 极小 P（含 P=0 下溢）原值保留、neg_log10_p=NA
6. project_summary 七类审计非占位；交付目录可独立移动
7. L1/L2/L3 回归不破
8. 全部发现与修复记录进 `docs/VALIDATION.md`；git 提交（信息格式：`fix/test: <描述>（真实数据验收）`）

## 8. 关键文件速查（调试入口）

| 关注点 | 文件 |
| --- | --- |
| 模块开关/目标计算 | `rules/utils/gwas_targets.py`（纯函数；改收集器必须同 PR 改 `rules/utils/datadeliver.py` 与对应规则 output，逐字一致） |
| 配置校验失败 | `schema/config.schema.yaml` + `rules/utils/validate.py` |
| 输入被拒 | `scripts/verify_inputs.py`（错误信息含具体文件与重复值行号） |
| 规范化失败 | `logs/02.input_gate/normalize.log`（含完整 bcftools 命令行） |
| QC 阻断 | `scripts/genotype_qc.py` fail() 输出的各阶段流失明细 |
| GEMMA 异常 | `07.association/<trait>/gemma.run.log` + `run_metadata.json`（status/exit_code/bindings） |
| 统计对账失败 | `scripts/harmonize_gemma.py` 结构错误输出（含冲突行明细，退出 2） |
| 交付校验 | `scripts/deliver_results.py`（MD5 不一致会列出期望/实际） |

## 9. 禁止事项（交接红线）

- 不得修改 `flow/FlowFrame/`（上游只读参考仓库）
- 不得删除或绕过 schema 的未实现功能拒绝（loco/farmcpu 等启用必须报错）
- 不得用 mock/合成数据结果充当 L4/L5 验收结论
- 不得在无真实软件版本台账时宣称"已验收"
- 不得提交分析产物（workflow/data_deliver 目录均被 .gitignore 排除；只提交代码/配置/文档/测试）
- 不得把 partial 写成 complete、不得事后缩小 family_size/M

## 10. 快速命令备忘

```bash
# dry-run（Step 0）
snakemake -s snakefile --config analysisyaml=<yaml> -n

# only_qc（Step 1）
snakemake -s snakefile --config analysisyaml=<yaml> --cores 8 --use-conda

# 单目标重跑（调试用）
snakemake -s snakefile --config analysisyaml=<yaml> 08.results/<trait>/gwas.tsv --rerun-triggers mtime

# 失败后恢复（不清空既有产物；Snakemake 默认行为 + 原子写入模式保证）
snakemake -s snakefile --config analysisyaml=<yaml> --cores 8 --rerun-incomplete

# 回归三件套
.testenv/bin/python -m pytest tests/test_unit.py -q
bash tests/test_dag.sh
bash tests/test_behavior.sh
```
