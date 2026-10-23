# GWASFlow 验证报告（VALIDATION.md）

版本：v0.1.0；日期：2026-10-23；验证环境：macOS 15 + Python 3.12.13 + Snakemake 9.27.1 + pytest 9.1.1（.testenv）。

## 1. 分层验证结果

| 层 | 手段 | 结果 | 说明 |
| --- | --- | --- | --- |
| L1 单元 | pytest（tests/test_unit.py） | **44/44 通过** | 覆盖验收场景 1/2/4/5/6/7/8/16/17 + CLI 横切契约（--help 快照、错误码 0/2/3、输出目录非空保护） |
| L2 DAG | snakemake dry-run（tests/test_dag.sh） | **21/21 通过** | 覆盖场景 3/10/12/14/15；基线 DAG 无循环依赖；两性状 prepare/run_gemma 作业数=2（无相互依赖边）；only_qc 截断五断言；loco=true Schema 层拒绝 |
| L3 行为 | mock 外部软件端到端（tests/test_behavior.sh） | **25/25 通过** | 全链 15 规则实际运行；场景 9 部分失败传导（failed_software/partial/健康性状不受影响）；场景 11 交付独立移动校验；失败恢复（QC 产物保留）；effective_config 覆盖生效；门控 baseline_only |
| L4 外部软件对账 | tests/commands/ 独立命令 | **未执行（M3）** | 真实 PLINK2/GEMMA/BCFtools 二进制与版本未定；记录为未验证项 |
| L5 真实数据 | — | **未执行（M4）** | 真实植物群体数据缺失；记录为未验证项 |

## 2. 关键行为验证明细（L3 端到端）

- **产物链齐备**：01.gate（三件+effective_config）→ 02.genotype（三件套+variant_map+审计）→ 03.genotype_qc（12 件）→ 04.structure → 05.phenotype（每性状 8 件）→ 06.kinship（K+顺序+验证）→ 07.association（每性状 5 件+MAC）→ 08.results（每性状五件套）→ 09.loci → 10.local_ld → 11.candidates → 12.plots → 13.summary → 14.data_deliver（三件套+copy 交付树）。
- **部分失败**：leaf_width 注入 GEMMA 非零退出 → 失败性状产出合法空表 + status=failed_software；leaf_length 正常完成；project_summary 整体状态 partial（不伪装 complete）。
- **交付独立移动**：交付目录移动后 delivery_manifest.md5 校验通过；analysis_results 树内容完整（copy 模式）。
- **失败恢复**：表型输入变化触发重算，既有 QC 产物（cohort.fam MD5）不变。
- **参数一致性**：traits.tsv 单点派生（n_pcs/协变量/transform）；analysisyaml 覆盖值反映进 effective_config.yaml。

## 3. 参数三方一致性核对（schema ↔ SCRIPT_PARAMETERS ↔ 命令构建器）

| 参数 | 单位 | schema | 手册 | 命令构建器 |
| --- | --- | --- | --- | --- |
| pre_geno/mind/geno/maf | 比例 0..1（maf≤0.5） | ratio01/maximum 0.5 | §3 比例 | genotype_qc.py --pre-geno 等 float |
| mac | 正整数 ≥1 | minimum 1 | §3 计数 | --mac int |
| prune_window_kb/clump_kb | kb 正数 | exclusiveMinimum 0 | §10 kb | --prune-window-kb/--clump-kb |
| prune_step | 标记数正整数 | minimum 1 | §10 标记数 | --prune-step int |
| flank_bp/fallback_flank_bp | bp 整数 ≥0 | minimum 0 | §10 bp | --flank-bp/--fallback-flank-bp int |
| chr_count | 染色体对数 1..95 | 1..95 | §3 1..95 | --chr-count int |
| pca_components/n_pcs | 维数 | minimum 1 / 非负 | §2/§3 | --pca-components/--n-pcs int |
| alpha | (0,1) Decimal 字符串 | pattern | §5 | --alpha Decimal |
| family_size | 正整数 | minimum 1 | §5 | --family-size int |
| primary_test | wald/lrt/score/all | enum | §4 | --test choices |
| kinship.type | centered/standardized | enum | §9.3 | --type choices |

各脚本 `--help` 输出与手册契约一致（L1 快照测试）。

## 4. 未验证项清单（不冒称已验证）

| 项 | 原因 | 后续动作 |
| --- | --- | --- |
| M1 参考包脚本重实现与参考行为差异 | 参考包 zip 缺失 | 以参数手册契约为唯一依据；L1 覆盖文档声明的 8 个合成测试场景 |
| M2 生菜组装真实版本/染色体条数 | 未提供 | 全配置化占位；部署前核对 |
| M3 真实软件数值对账（场景 13/18） | 二进制与版本未定 | 固定版本后执行 tests/commands/ 对照 |
| M4 真实群体端到端 | 数据缺失 | 上线前校准（含结构/亲缘模拟） |
| M5 rich-loguru 插件生产接线 | 本地缺失 | 三级降级链兜底；生产安装插件 |
| M6 BioReport 渲染器 | 本地不可用 | report=false；启用即显式报错 |
| M7 种子策略 | 未定 | seed 显式默认值进 project_summary |

## 5. 已知限制

- 统一统计内存排序（小/中规模）；大规模流式为后续扩展。
- define_loci 的 clump 成员边界解析依赖 plink2 .clumped 格式（解析失败回退固定窗口并记录 fallback_reason）。
- FAI offset 为占位（mock 场景；真实 bcftools 自行索引）。
- R 绘图在无 Rscript 环境降级为 mock（L3）；真实 R 验收待 M3。
