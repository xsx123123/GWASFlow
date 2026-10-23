# FLOWFRAME VARIANT PROFILE：input_profile=variants

版本：v0.1.0；日期：2026-10-23。
状态：**框架扩展提案**（在被上游 [CygnusX-Labs/FlowFrame](https://github.com/CygnusX-Labs/FlowFrame) 采纳前，不标注为"满足现行 FASTQ 模板全部条款"）。

## 1. 动机

FlowFrame v2.0 严格模板的 02/03 步（file_convert_md5 + short_read_qc）针对 **FASTQ + FastQC** 设计。GWASFlow 的输入是**群体基因型**（VCF/BCF/BED），不存在 R1/R2 概念，也不适用 FastQC 短读质检。

本声明定义 `input_profile=variants` 的文件级校验与基因型 QC 语义，替代 FASTQ 模板的对应步骤；**禁止**为 VCF 伪造 R1/R2 文件名或调用 FASTQ 预处理工具假装完成基因型标准化。

## 2. variants profile 的门控语义（对应 FASTQ 模板 02 步）

| FASTQ 模板条款 | variants profile 对应 |
| --- | --- |
| FASTQ 命名识别（R1/R2 软链） | 基因型格式识别（vcf/bcf/bed）与组件完整性（VCF/BCF 必含 .csi/.tbi 索引；BED 三件套齐备） |
| seq_preprocessor 建链 + MD5 | verify_inputs.py 逐文件角色/路径/可读性校验 + MD5 三态（verified/FAIL/baseline_only） |
| json_md5_verifier 比对 + TSV 报告 | md5_check.tsv（Status 列）+ input_gate.json（绑定配置 SHA-256 与输入 MD5） |
| raw_data_md5.json 审计层 | normalize_audit.json（多等位/INDEL 排除计数、FILTER 策略语义、DP/GQ 执行状态、染色体映射审计） |

门控接线七条军规原样继承：门控产物为硬门控（显式 input）；输入文件进 input 不进 params；规则 log 与工具自身 --log-file 分离。

## 3. variants profile 的 QC 语义（对应 FASTQ 模板 03 步）

- **分阶段基因型 QC**（四阶段独立执行，禁止单混合命令）：初步位点缺失 → 样本缺失 → 最终位点缺失/MAF/MAC → 最终队列重新统计。
- **非人类语义固定**：`--chr-set {N} no-x no-y no-xy no-mt --nonfounders`；HWE 默认关闭（自交/近交植物）；高亲缘材料保留。
- **基因型规范化前置**：REF 对 FASTA 严格校验（--check-ref e）、双等位 SNP 筛选、稳定 variant_id、染色体重命名（未知染色体报错，禁止压 0）。

## 4. 与框架目录规范的差异声明

| 目录 | 框架规范 | GWASFlow 现状 | 原因与降级方案 |
| --- | --- | --- | --- |
| `src/` | JZ_Tools submodule（日志/交付工具） | 暂不建立 | M5：本地未提供；降级为三级日志链（D9）+ Python 自研交付器（D6）；README 标注生产必须安装插件 |
| `report/` | bioreport submodule | 暂不建立 | M6：BioReport 本地不可用；报告能力为"结构化摘要 + 静态图表"；`report=true` 在渲染器验收前显式报错退出 2 |
| `rules/manifests/`、`agent/` | FlowAgent 入库资产 | 暂不建立 | FlowAgent 完整接入属后续工作包；不用空壳元数据冒充可用能力 |

## 5. 目标生成的 trait 维度扩展

FASTQ/RNA 模板以 sample×contrast 维度展开目标；GWASFlow 引入 **trait 维度**（`{trait}` wildcard，04/07/08/09/10/11/12 规则按性状独立展开，性状间无依赖边）。TargetPlan 纯函数契约（不改 config/无 IO/无 sleep/唯一调用点）与框架规范 §5.6 一致，并修复了模板 DataDeliver 的副作用缺陷（改写全局 config、time.sleep）。

## 6. 后续工作包边界（本轮不实现，启用必须明确报错）

LOCO（按染色体排除建 K）、FarmCPU/BLINK（GAPIT 独立进程）、基因型填充、多环境联合模型、人类 GWAS profile（HWE/性别/祖源/病例对照）、PGEN/HapMap 导入、流式大规模统计引擎（BGZF/Tabix/Parquet）、FlowAgent 完整入库（manifests/ports/catalog）、BioReport 完整 GWAS 渲染器。以上能力在 schema enum 与 `validate_module_switches` 双层拒绝（D8），错误信息指明"属后续工作包"。

## 7. 上游采纳路径

本 profile 作为独立扩展声明维护于 GWASFlow 仓库；上游 FlowFrame 若采纳 variants 输入类型，应：
1. 在框架规范的输入 profile 枚举中登记 `variants`；
2. 将本声明 §2/§3 的门控与 QC 语义合并入严格模板的 02/03 步通用条款；
3. 在 FlowAgent 入库手册补充群体基因型端口说明（组装/倍性/硬分型-剂量/样本顺序/性状变换/效应语义）。
