# GWASFlow

来源声明：本工程基于 [FlowFrame](https://github.com/CygnusX-Labs/FlowFrame)（commit `24a0cd3e`，框架规范 v2.0）构建，遵循 Flow 家族工程规范（入口九步初始化、配置四件套、TargetPlan 纯函数目标生成、资源档位、日志降级链与复制交付契约）。框架扩展声明见 [docs/FLOWFRAME_VARIANT_PROFILE.md](docs/FLOWFRAME_VARIANT_PROFILE.md)。

植物二倍体、硬分型双等位 SNP、连续数量性状的 GEMMA LMM 关联分析流程（多性状独立分析）。生菜（Lactuca sativa）为重点适配对象，染色体对数与参考版本全配置化。

## 适用边界

- **支持**：VCF/BCF（含索引）或 BED 三件套（+已核验 variant_map）、材料级连续数量性状、Wald/LRT/score 预声明主检验、全队列 K + 性状缺失掩码路线。
- **不支持（启用即报错，属后续工作包）**：LOCO、FarmCPU/BLINK、基因型填充、多环境联合模型、人类 GWAS profile、PGEN/HapMap 导入、大规模流式统计引擎、完整 FlowAgent 入库、BioReport 渲染器。
- **规模边界**：统一统计首版面向小/中规模（内存排序）；大规模流式（BGZF/Tabix/Parquet）为后续扩展。

## 目录结构

```text
GWASFlow/
├── snakefile                 # 入口（小写）：九步初始化
├── config/                   # 配置四件套（config/reference/run_parameter/cluster_config）
├── schema/config.schema.yaml # 全键校验 + 未实现功能拒绝分支（D8 第一层）
├── rules/                    # 01..15 模块化规则（编号≠执行顺序）
│   └── utils/                # TargetPlan 纯函数、收集器、校验工具
├── scripts/                  # 14 个 Python CLI + 2 个 R CLI（参数手册契约）
├── envs/                     # conda 环境（版本占位待 M3 固定）
├── examples/                 # 配置覆盖示例 + 七类表格模板
├── tests/                    # L1/L2/L3 分层测试（覆盖 18 验收场景）
└── docs/                     # 实现文档/参数手册/框架对齐/variants profile
```

## 快速开始

### 1. 环境准备

```bash
# 测试/运行虚拟环境
uv venv .testenv --python 3.12
uv pip install --python .testenv/bin/python pytest pyyaml "snakemake>=9.9.0"
# 外部软件（M3：固定版本后重验收）
conda env create -f envs/gwasflow-genotype.yaml   # python+bcftools+plink2
conda env create -f envs/gemma.yaml               # python+gemma+numpy
conda env create -f envs/r-plot.yaml              # R+ggplot2
```

### 2. 配置

复制 `examples/analysisyaml.example.yaml` 为项目配置，替换绝对路径与 Genome_Version：

```bash
cp examples/analysisyaml.example.yaml /data/project/01.workflow/analysis.yaml
# 编辑：输入表路径、基因型输入、workflow/data_deliver、参考资源路径
```

### 3. 框架初始化验证（dry-run）

```bash
snakemake -s snakefile --config analysisyaml=/data/project/01.workflow/analysis.yaml -n
# 生效模块开关清单与自动启用原因会打印到控制台/日志
```

### 4. 实际运行

```bash
# 本地
snakemake -s snakefile --config analysisyaml=/data/project/01.workflow/analysis.yaml \
    --cores 8 --use-conda --rerun-triggers mtime

# 集群（示例 slurm）
snakemake -s snakefile --config analysisyaml=... --executor slurm --use-conda \
    --default-resources slurm_account=your_account
```

### 5. 模块级目标请求

```bash
# 仅重跑某性状统一统计
snakemake -s snakefile --config analysisyaml=... 08.results/leaf_length/gwas.tsv
# only_qc 模式（analysisyaml 设置 only_qc: true）
```

## 测试

```bash
.testenv/bin/python -m pytest tests/test_unit.py -v   # L1 单元（44 项）
bash tests/test_dag.sh                                 # L2 DAG（21 项）
bash tests/test_behavior.sh                            # L3 mock 端到端（25 项）
```

覆盖矩阵见 `tests/coverage_matrix.md`（18 个验收场景无遗漏映射）。

## 重要声明

- **QC 默认阈值（mind=0.10/geno=0.10/maf=0.05/mac=5）为审阅起点，非普适植物标准**；自交/近交植物默认关闭 HWE 过滤、不自动删除高亲缘材料。
- **生产环境必须安装 `snakemake-logger-plugin-rich-loguru`**（M5：本地缺失，代码含三级降级链兜底，运行时附 `--logger rich-loguru`）。
- 生菜示例值（染色体 9 对、`Lsat_Salinas_v11`）为占位（M2）：部署前必须核对真实组装。
- 外部软件数值验收（GEMMA 对账）记录为未验证项（M3），见 `docs/FRAMEWORK_ALIGNMENT.md` §6。

## 文档索引

| 文档 | 内容 |
| --- | --- |
| [docs/GWASFlow_IMPLEMENTATION.md](docs/GWASFlow_IMPLEMENTATION.md) | 核心实现文档（链路/规则清单/方法约束） |
| [docs/SCRIPT_PARAMETERS.md](docs/SCRIPT_PARAMETERS.md) | 脚本参数与软件命令手册（CLI 契约） |
| [docs/FRAMEWORK_ALIGNMENT.md](docs/FRAMEWORK_ALIGNMENT.md) | 框架对齐调查与 M1–M7 缺失信息 |
| [docs/FLOWFRAME_VARIANT_PROFILE.md](docs/FLOWFRAME_VARIANT_PROFILE.md) | input_profile=variants 框架扩展声明 |
| tests/README.md | 分层测试运行入口与未验证项清单 |
