# GWASFlow 框架对齐调查报告（FRAMEWORK_ALIGNMENT.md）

版本：v0.1.0-survey；日期：2026-10-23；性质：需求规格阶段的调查结论，非设计文档。
配套需求规格：`.codeartsdoer/specs/gwas_gemma_loop/spec.md`。

## 1. 参考仓库实际提交记录

| 资源 | 位置 | 实际提交 / 状态 |
| --- | --- | --- |
| FlowFrame | 本地 `flow/FlowFrame` | `24a0cd3e8be92b744a4c6f9df7929298d2318b0c`（git HEAD，与实现文档 §3 记录一致） |
| FlowFrame 规范 | `flow/FlowFrame/flow-framework/references/flow_framework_spec.md` | v2.0，1683 行 |
| FlowFrame 技能 | `flow/FlowFrame/flow-framework/SKILL.md` | 含 RNAFlow 提取的全部模板（`assets/templates/`） |
| FlowAgent | `flow/FlowFrame/FlowAgent/` | TypeScript 工程（src/docs/prompts/eval），入库手册在 `docs/planning/` |
| RNAFlow | 本地无副本 | 仅 GitHub 引用；通用层模板已内置于 FlowFrame `assets/templates/`，可作生产参考 |
| GWASFlow 参考包 | `GWASFlow_Implementation_Package_20261023.zip` | **未在本地找到**（见 §3） |

## 2. 当前仓库状态

- 工作目录 `/Users/zj/zj_code/GWASFlow`：仅有 `docs/`（两份核心文档）、`flow/FlowFrame`（上游参考框架）、`.codeartsdoer/`。
- **无既有 GWASFlow 实现代码**，本轮按全新工程创建；不存在需要保留的用户改动。
- 上游 FlowFrame 仅作参考，框架扩展写入本项目文档，不修改上游仓库。

## 3. 参考包可用性与复用结论

参考包 zip 未找到，本地仅存其两份文档：

| 参考包内容 | 文档声称状态 | 本地实际 | 复用结论 |
| --- | --- | --- | --- |
| `docs/GWASFlow_IMPLEMENTATION.md` | 已提供 | ✅ 存在（355 行） | 作为需求与验收依据直接采用 |
| `docs/SCRIPT_PARAMETERS.md` | 已提供 | ✅ 存在（502 行） | 作为参数契约（CLI、错误码、单位）直接采用 |
| `scripts/prepare_trait_inputs.py` | 已实现+合成测试 | ❌ 缺失 | 按参数手册 §2 契约重新实现 |
| `scripts/harmonize_gemma.py` | 已实现+合成测试 | ❌ 缺失 | 按参数手册 §5 契约重新实现 |
| `scripts/annotate_loci_genes.py` | 已实现+合成测试 | ❌ 缺失 | 按参数手册 §6 契约重新实现 |
| `commands/01_genotype_qc.sh` | 语法通过 | ❌ 缺失 | 按参数手册 §3 契约重新实现（生产应拆分为可恢复规则） |
| `commands/02_run_gemma_reference.sh` | 语法通过 | ❌ 缺失 | 按参数手册 §4 契约重新实现（仅作独立对照） |
| tests（8 个合成测试） | 已通过 | ❌ 缺失 | 需求中重新声明测试场景（spec §5 各模块验收条件） |
| examples 表格 | 已提供 | ❌ 缺失 | 按实现文档 §4 契约重建 |

**结论**：参考代码不能直接复用；所有脚本按文档契约重新实现并通过文档声明的合成测试验收。参考代码本身即使可得也需审查而非直接投产。

## 4. FlowFrame 对齐要点（需求层面）

本轮 GWASFlow 必须对齐的框架能力（详细业务规则见 spec.md §5.12）：

1. **入口**：小写 `snakefile` 8 步初始化、四件套 configfile、`analysisyaml` 最高优先级覆盖、Schema 强制校验、workdir 重定向。
2. **规则分层**：通用模板层（01.common → 02 输入门控 → 03 输入 QC）→ 组学专用层（04..N-2）→ 收尾层（deliver/Report）。
3. **门控**：所有直接消费外部输入的规则显式挂门控产物于 `input`；输入文件必须进 `input` 不进 `params`。
4. **目标生成**：不可变 TargetPlan 纯函数一次计算；AnalysisTargets/DataDeliver 复用同一计划；deliver/report 产物必须进入 rule all；模块开关显式默认值，开关静默失效按 P0 处理。
5. **资源与日志**：rule_resource 档位、conda envs 固定版本、日志编号目录且不与工具自身日志双写。
6. **交付**：默认 copy 模式，manifest 记录文件清单与校验和，交付目录可独立移动。

## 5. 框架扩展声明（需求来源）

- 现行严格模板 02/03 步仅针对 FASTQ + FastQC。GWASFlow 输入为群体基因型（VCF/BCF/BED），必须定义 `input_profile=variants` 的文件级校验与基因型 QC（对应 `docs/FLOWFRAME_VARIANT_PROFILE.md`，由设计/实施阶段产出）。
- 该扩展在被上游采纳前**不标注为“满足现行 FASTQ 模板全部条款”**；禁止为 VCF 伪造 R1/R2 文件名或调用 FASTQ 预处理工具。
- 框架中旧地址 `xsx123123/FlowFrame` 与现行 `CygnusX-Labs/FlowFrame` 并存，README 来源声明按现行框架处理并统一迁移。

## 6. 明确缺失信息（需用户/部署方提供）

| 编号 | 缺失项 | 影响 | 当前处理 |
| --- | --- | --- | --- |
| M1 | 参考包 zip 内脚本/测试源码 | 无法直接复用，只能按文档契约重写 | 记录为重实现（§3） |
| M2 | 生菜参考组装真实版本、染色体命名与条数（文档示例 `Lsat_Salinas_v11`/9 对） | 染色体映射与注释边界无法实证 | 染色体编码、参考版本全部配置化，示例值仅占位 |
| M3 | PLINK 2 / GEMMA / BCFtools 二进制与精确版本 | 软件过滤对账、数值一致性验收无法执行 | 验收分层：可独立实施部分先完成，外部软件验收记录为未验证项 |
| M4 | 真实生菜（或等价植物）群体数据 | 端到端真实运行验收无法执行 | 小型合成数据验证 DAG 与行为，真实数据验收记录为未验证项 |
| M5 | JZ_Tools 日志插件子模块（`src/`）本地未提供 | `--logger rich-loguru` 强制项无法接线 | 设计阶段决定降级方案（标准 Snakemake 日志），记录为待补 |
| M6 | BioReport 仓库本地不可用 | GWAS 报告适配无法确认接入条件 | 报告能力声明为“结构化摘要+静态图表”；BioReport 适配为条件性最小实现 |
| M7 | 合成数据与真实分析的种子策略 | 复现记录中种子字段 | 配置化 seed 字段，默认值显式声明 |

## 7. 后续阶段交接

- 本文档与 `.codeartsdoer/specs/gwas_gemma_loop/spec.md`（需求规格）共同构成下一阶段（spec-design-agent 实现方案创建）的输入。
- 本轮范围：植物二倍体、硬分型双等位 SNP、连续数量性状的 GEMMA 主分析闭环（多性状独立分析）。LOCO、FarmCPU/BLINK、填充、多环境、人类 GWAS、FlowAgent 完整接入为独立工作包，启用时必须明确报错。
