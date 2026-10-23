# GWASFlow

来源声明：本工程基于 [FlowFrame](https://github.com/CygnusX-Labs/FlowFrame)（commit `24a0cd3e`，框架规范 v2.0）构建，遵循 Flow 家族工程规范（入口九步初始化、配置四件套、TargetPlan 纯函数目标生成、资源档位、日志降级链与复制交付契约）。框架扩展声明见 [docs/FLOWFRAME_VARIANT_PROFILE.md](docs/FLOWFRAME_VARIANT_PROFILE.md)。

植物二倍体、硬分型双等位 SNP、连续数量性状的 GEMMA LMM 关联分析流程（多性状独立分析）。染色体对数、参考版本与模块开关全配置化。

**当前状态（v0.1.0）**：已在拟南芥 AtPolyDB 真实数据（194 样本 × 214,051 SNP，FT10 性状）上完成端到端验收——全链运行、独立 GEMMA 数值对账（口径 A 逐位一致）、QC 剔除集对账（口径 C 一致）、FT10 已知主效位点（FRI 区域）生物学复现，且 `--use-conda` 与本机 PATH 两种部署模式结果一致。详见 [docs/VALIDATION.md](docs/VALIDATION.md) 与 [docs/ATPOLYDB_TEST_GUIDE.md](docs/ATPOLYDB_TEST_GUIDE.md)。

## 适用边界

- **支持**：VCF/BCF（含索引）或 BED 三件套（+已核验 variant_map）、材料级连续数量性状、Wald/LRT/score 预声明主检验、全队列 K + 性状缺失掩码路线、静态图（svg/pdf/png）与交互图（plotly WebGL HTML）。
- **不支持（启用即报错，属后续工作包）**：LOCO、FarmCPU/BLINK、基因型填充、多环境联合模型、人类 GWAS profile、PGEN/HapMap 导入、大规模流式统计引擎、完整 FlowAgent 入库、BioReport 渲染器。
- **规模边界**：统一统计首版面向小/中规模（内存排序，实测 214k 位点正常）；大规模流式（BGZF/Tabix/Parquet）为后续扩展。

## 目录结构

```text
GWASFlow/
├── snakefile                 # 入口（小写，可自动发现）：九步初始化
├── config/                   # 配置四件套（config/reference/run_parameter/cluster_config）
├── schema/config.schema.yaml # 全键校验 + 未实现功能拒绝分支（D8 第一层）
├── rules/                    # 01..15 模块化规则（编号≠执行顺序；docstring 均为 Input/Process/Output 三段式）
│   └── utils/                # TargetPlan 纯函数、收集器、校验工具
├── scripts/                  # 16 个 Python CLI + 2 个 R CLI（参数手册契约）
├── envs/                     # conda 环境声明（7 个；--use-conda 时按规则自动创建）
├── examples/                 # 配置覆盖示例 + 七类表格模板
├── tests/                    # L1/L2/L3 分层测试（覆盖 18 验收场景）
└── docs/                     # 实现文档/参数手册/框架对齐/测试指南/验收记录
```

## 快速开始

### 1. 环境准备

**调度环境**（snakemake + 测试）：

```bash
uv venv .testenv --python 3.12
uv pip install --python .testenv/bin/python pytest pyyaml "snakemake>=9.9.0"
```

**外部工具**，两种方式二选一：

- **方式 A（推荐）：`--use-conda` 自动创建**——每条规则声明了 `conda:` 指令，运行时自动建 5 个环境（gwasflow-python / plink2 / gemma / r-plot / plotly，约 2.2GB）。前置依赖：

  ```bash
  conda install -n base -c conda-forge mamba        # --conda-frontend mamba 需要
  uv pip install --python .testenv/bin/python \
      /path/to/snakemake-logger-plugin-rich-loguru  # --logger rich-loguru 需要（本地源码或 PyPI）
  ```

- **方式 B：本机 PATH 直连**——plink 1.9、plink2、gemma、bcftools、R（ggplot2 + jsonlite + svglite）、Python 3.10+ 装入 PATH 即可，运行时不加 `--use-conda`。

两种模式结果一致（已实测：gwas.tsv 指纹相同）。

### 2. 配置（analysisyaml）

复制 `examples/analysisyaml.example.yaml` 为项目配置，放进项目工作目录（如 `/data/project/01.workflow/analysis_config.yaml`），编辑：输入表路径、基因型输入（`genotype_bfile`/`genotype_input`）、`workflow`/`data_deliver`、`GWAS_reference` 资源路径。绝对路径传入，位置不限。真实数据模板可参考 `tests/atpolydb/analysisyaml.atpolydb_smoke.yaml`。

### 3. dry-run 验证

```bash
snakemake --cores 100 --conda-frontend mamba --use-conda --logger rich-loguru \
    --config analysisyaml=/data/project/01.workflow/analysis_config.yaml -n
# 生效模块开关清单与自动启用原因会打印到控制台/日志
```

### 4. 实际运行

```bash
# 本地（多核，推荐提交格式）
snakemake --cores 100 --conda-frontend mamba --use-conda --rerun-triggers mtime \
    --logger rich-loguru \
    --config analysisyaml=/data/project/01.workflow/analysis_config.yaml \
    --rerun-incomplete

# 集群（示例 slurm；集群执行档位见 config/cluster_config.yaml）
snakemake --cores 100 --executor slurm --use-conda --conda-frontend mamba \
    --logger rich-loguru \
    --config analysisyaml=/data/project/01.workflow/analysis_config.yaml \
    --default-resources slurm_account=your_account
```

说明：`snakefile` 为小写且可自动发现，无需 `-s`；`--logger rich-loguru` 需装插件（未装时代码有三级降级链兜底，仅提示不阻断）；conda 环境按项目缓存在 `<workflow>/.snakemake/conda/`。

### 5. 模块级目标请求

```bash
# 仅重跑某性状统一统计
snakemake --cores 8 --config analysisyaml=... 08.results/leaf_length/gwas.tsv
# only_qc 模式（analysisyaml 设置 only_qc: true，只跑门控+规范化+QC）
```

### 产物速览

全链产出 14 个目录：QC 审计（`03.genotype_qc/`）、PCA/K（`04`–`06`）、GEMMA 关联与统一统计（`07`–`08`）、位点与区域 LD（`09`–`10`）、静态图 + **交互式曼哈顿/QQ 图**（`12.plots/{trait}/manhattan_interactive.html`，plotly WebGL 全量渲染，离线可开）、项目摘要与交付树（`13`–`14`）。

## 测试

```bash
.testenv/bin/python -m pytest tests/test_unit.py -q    # L1 单元（46 项）
bash tests/test_dag.sh                                 # L2 DAG（21 项）
bash tests/test_behavior.sh                            # L3 mock 端到端（25 项，含交互图规则）
```

覆盖矩阵见 `tests/coverage_matrix.md`（18 个验收场景无遗漏映射）；AtPolyDB 数据包合同测试见 `docs/ATPOLYDB_TEST_GUIDE.md` Step 7。

## 验证状态（本机实测台账）

| 项 | 结果 |
| --- | --- |
| AtPolyDB 真实数据全链（smoke 10,703 SNP / 完整集 214,051 SNP） | ✅ 均 100% 完成；FT10 在 FRI 区域复现 Bonferroni 显著主效位点（chr4:6,682,595，p=1.15e-7） |
| L4 数值对账（口径 A：独立 GEMMA 同参数复现） | ✅ smoke 逐位一致（偏差 0）；完整集 203,813 值中 1 个 P 值相对偏差 4.77e-07（GEMMA 运行间浮点抖动，<1e-6 容差） |
| QC 对账（口径 C：四阶段剔除集） | ✅ 逐阶段集合一致 |
| `--use-conda` vs PATH 部署模式 | ✅ gwas.tsv 指纹一致 |
| 软件版本（实测） | GEMMA 0.98.5、PLINK 1.9.0-b.7.11、PLINK 2.0.0-a.7.x、R 4.5.x、Python 3.12/3.14 |

详细对账记录与未验证项收缩见 [docs/VALIDATION.md](docs/VALIDATION.md)。

## 重要声明

- **QC 默认阈值（mind=0.10/geno=0.10/maf=0.05/mac=5）为审阅起点，非普适植物标准**；自交/近交植物默认关闭 HWE 过滤、不自动删除高亲缘材料。
- 生菜示例值（染色体 9 对、`Lsat_Salinas_v11`）为占位（M2）：部署前必须核对真实组装。
- BED 路线 variant_map 契约：`variant_id` 必须登记为 BIM 的 SNP ID（外部原始 ID）；BIM A1/A2 为观测等位基因，无 REF/ALT 断言。
- 交互图为探索用途（HTML 全量渲染，超大数据集可配 `plotting.max_render_points` 抽样）；静态三格式图为发表用途。
- `15.Report` 渲染器（M6）仍未验收：`report: true` 会显式报错退出 2，当前交付物为结构化摘要 + 静态/交互图表。

## 文档索引

| 文档 | 内容 |
| --- | --- |
| [docs/GWASFlow_IMPLEMENTATION.md](docs/GWASFlow_IMPLEMENTATION.md) | 核心实现文档（链路/规则清单/方法约束） |
| [docs/SCRIPT_PARAMETERS.md](docs/SCRIPT_PARAMETERS.md) | 脚本参数与软件命令手册（CLI 契约） |
| [docs/FRAMEWORK_ALIGNMENT.md](docs/FRAMEWORK_ALIGNMENT.md) | 框架对齐调查与 M1–M7 缺失信息 |
| [docs/FLOWFRAME_VARIANT_PROFILE.md](docs/FLOWFRAME_VARIANT_PROFILE.md) | input_profile=variants 框架扩展声明 |
| [docs/ATPOLYDB_TEST_GUIDE.md](docs/ATPOLYDB_TEST_GUIDE.md) | AtPolyDB 真实数据 L4/L5 验收指南（分步命令） |
| [docs/HANDOVER.md](docs/HANDOVER.md) | 通用交接文档（L1–L5 验收路径与风险区） |
| [docs/VALIDATION.md](docs/VALIDATION.md) | 验收记录（对账结果/版本台账/未验证项） |
| tests/README.md | 分层测试运行入口与未验证项清单 |
