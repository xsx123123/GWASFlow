# GWASFlow 详细实现文档

版本：v0.1.0-reference  
日期：2026-10-23  
目标：为 Flow 家族新增植物群体 GWAS 工作流，供开发、方法审阅和验收使用。

## 1. 交付范围与实现状态

本文将“目标实现”与“本参考包已实现”分开标注。当前包包含工程规范、脚本参数、配置模板及参考工具，不是已经完成的生产 GWASFlow。

| 内容 | 当前状态 |
| --- | --- |
| 架构、输入输出、规则接线、开发与验收规范 | 已提供 |
| prepare_trait_inputs.py | 已实现并通过合成数据测试 |
| harmonize_gemma.py | 已实现并通过合成数据测试 |
| annotate_loci_genes.py | 已实现并通过合成数据测试 |
| PLINK QC/PCA、GEMMA K/LMM 命令封装 | 已提供；Shell 语法通过，外部软件未执行 |
| 配置四件套及样本/性状等表格 | 已提供模板；路径为占位值 |
| 生产 snakefile、全部 Snakemake rules、Schema | 待按本文实现 |
| VCF 输入门控与规范化、完整绘图/报告 | 待实现 |
| LOCO、FarmCPU/BLINK、FlowAgent 入库资产 | 后续阶段 |
| 真实植物群体、集群和端到端验证 | 尚未完成 |

Python 脚本依赖 Python >=3.10 标准库，无 pandas/scipy 依赖。当前验证环境 Python 3.14.6。软件参数依据官方文档和源码核对，真正部署必须固定精确版本并重新验收。

## 2. 适用场景与输入边界

首版：二倍体植物自然群体/关联群体，硬分型双等位 SNP，连续数量性状。重点适配生菜材料，统计参数仍需结合真实群体审阅。

输入必须是同组装下联合分型或等价一致处理的群体基因型。不要把不同个体、不同参考组装的单样本 VCF 直接拼接成可分析群体。

支持目标：VCF/BCF、BED/BIM/FAM。PGEN 和 HapMap 作为独立导入适配，未验证前禁止宣称支持。剂量输入不能静默转为 BED；多倍体、二分类/等级/生存性状、多环境联合模型要单独建 profile。

FASTQ→比对→联合分型在上游 WGSFlow 等流程完成；本流程以可靠的群体基因型和表型为起点。

## 3. FlowFrame 对齐与需要扩展的部分

基线来源为 CygnusX-Labs/FlowFrame 的开发技能、完整规范及入库手册，并参考 RNAFlow 的公开入口和文档。此前浏览的 FlowFrame main 为 24a0cd3e8be92b744a4c6f9df7929298d2318b0c；正式建仓须记录重新取得的实际提交。

继承框架的入口、配置、目标生成、资源、日志和交付接口。GWAS 新增的计算维度是 cohort/trait/model，不能沿用“差异比较组=性状”的隐含假设。

当前严格模板的 02/03 步针对 FASTQ 和 FastQC。因此先写 docs/FLOWFRAME_VARIANT_PROFILE.md，声明 input_profile=variants，定义文件级校验和基因型 QC。它是框架扩展提案，在被采纳前不标注为“满足现行 FASTQ 模板全部条款”。

不要为 VCF 伪造 R1/R2 文件名，也不要调用 FASTQ 预处理工具假装完成基因型标准化。

生产仓库继续采用：snakefile、config/ 四件套、schema/config.schema.yaml、rules/utils/、scripts/、envs/、src/、report/、docs/、examples/。新增 rules/manifests/ 与 agent/，在真正接入 FlowAgent 时完成。

README 来源声明按现行框架处理，同时补充当前组织地址；框架技能中旧的 xsx123123/FlowFrame 地址需要统一迁移，不能把旧地址当成仍可访问的事实。

## 4. 数据身份与输入契约

### 4.1 样本身份

samples.csv 至少含 sample,fid,iid。sample 是材料内部唯一键；FID/IID 精确对应 PLINK 文件，二元组必须唯一。不要靠下划线拼接身份后再拆分。

可选 sample_name/group/include：显示名、群体元数据和显式是否纳入。group 不默认进入模型，不生成 Treat_vs_Control。phenotypes.tsv 和 covariates.tsv 使用 sample 关联。

VCF 导入的示例采用 --double-id，因此 FID=IID=VCF sample ID；本包 examples/samples.csv 对应这种编码。如果导入策略改为 FID=0，必须相应修改样本表。

### 4.2 表型和协变量

phenotypes.tsv：每材料一行，每性状一列，缺失为 NA/NaN/./空值。-9 不被 Python 准备脚本默认视为缺失，防止误删合法负值；需要该哨兵时显式加到 --missing-values。

traits.tsv：trait、phenotype_column、trait_type、unit、transform、primary_model、numeric_covariates、categorical_covariates、n_pcs。脚本参数由这张表解析生成，不重复在多个位置维护同一性状定义。

只接受材料级性状。不同年份/地点/重复测量需先采用明确的实验设计模型得到材料级观测或估计，并保留模型/精度信息。BLUE/BLUP 或调整均值不能不加区分地当成等精度独立观测。

数值协变量与分类协变量必须显式区分；参考脚本按纳入材料的分类水平字典序选择参考组并登记。研究中若需固定参考组，生产实现应增加声明字段而非依赖隐式排序。

### 4.3 参考与变异身份

reference.yaml 登记当前组装所需的 FASTA/FAI、染色体表、重命名映射、GFF3 和可选功能注释。仅校验已启用模块需要的当前组装资源。

variant_map.tsv 是规范化结果的身份表：

| 列 | 约束 |
| --- | --- |
| variant_id | 唯一稳定 ID；建议 original_chrom:pos:ref:alt |
| original_id | 供应方原始 ID，可为 . |
| chrom | 最终展示的原始参考染色体名 |
| engine_chrom | PLINK/GEMMA 内部数值编码 |
| pos | 1-based 正整数 |
| ref/alt | 经 FASTA 核验的 A/C/G/T 双等位 SNP |

BED 的 A1/A2 不能独立证明 REF/ALT。本包 harmonizer 要求外部提供已核验的 variant_map，不从 BIM 第六列盲猜参考碱基。

## 5. 推荐数据流

~~~mermaid
flowchart TD
  G["群体基因型"] --> Gate["文件与参考门控"]
  T["样本、表型、协变量"] --> Gate
  Gate --> QC["基因型质控"]
  Gate --> P["每性状有效样本与设计矩阵"]
  QC --> S["LD pruning 与 PCA"]
  S --> P
  QC --> K["全队列 kinship"]
  QC --> A["GEMMA 主分析"]
  P --> A
  K --> A
  A --> H["统一关联统计"]
  H --> L["显著区域、LD 与候选证据"]
  L --> E["摘要、报告和复制交付"]
~~~

### 5.1 全队列与性状集合

先生成最终 QC cohort，再由 prepare_trait_inputs 为每个性状产生完整 N 行的表型/协变量文件。无效材料保留行位置，表型为 NA；有效材料进入 keep.tsv。

参考命令以全 1 的 kinship_all.pheno.txt 强制 kinship 计算覆盖全队列，不读取 FAM 原有缺失表型。关联时用相同 N 顺序的 K 和带 NA 的性状向量，让软件提取有效个体。

“全队列估计 K 后取子集”与“每性状重新估计 K”是不同方法，本参考包采用前者。生产流程可缓存共享 K，参考 Shell 为独立验证方便每次重新计算，尚未优化共享缓存。

## 6. 规则级实现清单

| 文件 | 关键规则/函数 | 输入 | 输出 |
| --- | --- | --- | --- |
| 01.common.smk | 公共初始化 | 生效配置 | 日志上下文 |
| 02.input_gate.smk | verify_inputs、normalize_genotype | 基因型、表格、参考、提供方校验和 | input_gate.json、md5_check.tsv、规范化 VCF、variant_map |
| 03.genotype_qc.smk | genotype_qc | 规范化基因型和 gate | cohort.bed/bim/fam、统计、剔除原因 |
| 04.phenotype_qc.smk | prepare_trait_inputs | FAM、表型、协变量、PCs | 性状向量、设计矩阵、keep、metadata |
| 05.population_structure.smk | prune、pca | 最终 cohort | prune.in、eigenvec、eigenval |
| 06.kinship.smk | compute_kinship、validate_kinship | cohort、all-ones phenotype | K、样本顺序、验证结果 |
| 07.association.smk | run_gemma | cohort、K、性状向量、设计矩阵 | 原始 assoc、日志、run_metadata |
| 08.result_harmonization.smk | harmonize_gemma | assoc、variant_map、两个 metadata | gwas.tsv、显著/异常/未报告表 |
| 09.locus_definition.smk | define_loci | 主检验、阈值、有效材料基因型 | lead、loci、区域来源信息 |
| 10.local_ld.smk | regional_ld | 区域及对应材料基因型 | LD、区域图索引、状态 |
| 11.candidate_genes.smk | annotate_loci_genes | loci、GFF、可选功能表 | candidates.tsv 与来源 |
| 12.visualization.smk | plot_gwas、plot_diagnostics | 统一统计、PC、K、表型 | 图、绘图数据、图表清单 |
| 13.export_summary.smk | export_project_summary | 真实模块产物与状态 | project_summary.json、artifact_index.tsv |
| 14.deliver.smk | deliver_results | 明确的分析/报告文件集合 | 复制结果、manifest、校验和 |
| 15.Report.smk | render_gwas_report | project_summary 和索引 | HTML/PDF；模块不可用时明确失败或摘要降级 |

编号不是执行顺序：04 的某些规则可以依赖 05 产物。所有直接消费外部输入的规则显式挂 gate，文件放 input，不藏在 params。

## 7. 目标生成、模块开关与环路避免

实现不可变 TargetPlan，字段至少为 analysis_outputs、summary_outputs、report_outputs、delivery_outputs、all_outputs。

目标函数只依据已解析配置、样本和性状定义计算路径，不改全局 config，不扫描分析结果，不启动软件，不等待文件。AnalysisTargets/DataDeliver 保留兼容接口，但复用同一份已计算计划。

模块依赖建议：

| 模块 | 必要上游 |
| --- | --- |
| population_structure | genotype_qc |
| phenotype_qc（使用 PCs） | population_structure |
| kinship_analysis | genotype_qc |
| association_analysis | genotype_qc、phenotype_qc、kinship_analysis |
| locus_analysis | association_analysis |
| local_ld | locus_analysis、相应样本/标记集 |
| candidate_genes | locus_analysis、基因位置注释 |
| visualization | 所画图对应的真实产物 |

默认值在 config/config.yaml 显式给出，解析依赖后输出生效开关及自动启用原因。only_qc=true 时仅保留输入、基因型/表型质控及其摘要，不自动运行 PCA/K/关联。

报告依赖分析产物，最终交付可依赖报告，但不能把包含报告/交付自身的 ALL_TARGETS 再作为其输入。推荐分析→摘要→报告→最终交付；独立的 deliver/report 开关不能互相隐式关闭。

locus 数量运行后才知道：首版每性状批量产生区域索引和目录，目录中每文件都进入 manifest。需要逐 locus 并行时使用 Snakemake checkpoint，不能在纯目标函数中读取尚不存在的统计表。

## 8. VCF/BCF 预处理实施要求

1. 验证提供方文件校验和，有提供方基线才称 verified；无基线时计算并标 baseline_only。
2. 校验索引和文件内容关系，REF 对 FASTA，一致性失败时阻断。不能用自动换 REF/ALT 代替组装确认。
3. 首版保留原始双等位 SNP，排除多等位和 INDEL 并记数。若未来拆分多等位，须重新定义其他 ALT、GT/DS 的语义。
4. 明确 FILTER 策略：PASS-only 或 PASS+“.”，后者的“.”代表未评估，不能称已通过原始质量过滤。
5. DP/GQ 若启用先按单样本基因型置缺失，再计算缺失率；字段缺失时阻断该过滤或明确记录未执行。
6. 保留 original_id 后生成稳定 variant_id，生成 variant_map，再把染色体映射为 engine_chrom。
7. 转 BED 前确认只含二倍体硬分型；缺失/半缺失策略明确。VCF dosage 不能静默损失。
8. 采用固定、独立 QC 阶段：初步位点缺失→样本缺失→最终位点缺失/MAF/MAC→重新统计。不要用一条混合命令假定过滤会迭代。

具体 BCFtools/PLINK 参数及待开发 wrapper 接口见参数手册。包内 01_genotype_qc.sh 从已经规范化的 BED 开始，不替代本节预处理和门控。

## 9. 质控与方法约束

示例参数 mind=0.10、geno=0.10、maf=0.05、mac=5，都是审阅起点，不是普适植物标准。

自交/近交植物不默认 HWE 过滤，不默认删除高亲缘材料；杂合度、群体离群先标记并展示。人类的 HWE、性别、祖源/病例规则需要独立 profile。

PCA 只使用经 QC 的 LD-pruned 集合；关联分析使用完整最终标记集。PC 数量按研究计划声明，不能按哪个维度获得更多显著位点来调参。

每个性状有效材料改变后，MAF/MAC 和可测试性可能改变。GEMMA 会在分析材料上重新应用自身过滤，因此要核对实际 M。参考 Shell 不额外实施 trait-specific MAC；生产实现需由该性状材料计算精确 allele counts 并生成允许测试标记集。不能把队列级 MAC 当成所有性状的精确 MAC。

表型不自动补均值，不自动正态化，不自动删离群材料。参考脚本只支持 none/log/log1p，变换会改变效应单位。秩不足或残余自由度不足即阻断，避免不合法模型继续运行。

## 10. GEMMA 主模型、LOCO 与补充模型

首版主模型：

y = intercept + SNP effect + declared covariates + selected PCs + u + error

其中 u 的协方差由 K 指定。主检验预声明，Wald 为本包主路径；LRT/score/all 需对选定软件版本的输出另验。

K 类型 centered/standardized 必须记录，不能跨项目缓存混用。当前参考封装只做 genome_wide K；LOCO 后续单独验收。保守的实现方式是先生成“排除受检染色体”的建 K 标记集，再按染色体关联并统一归档。不要只给全局 K 的结果改标签。

GEMMA 源码确认：同时提供 -bfile 和 -p 时使用外部表型；因此本包保留 FAM 的身份顺序，用 -p 指向完整 N 行的 trait phenotype。-c 文件显式含截距。GEMMA 自身的 -maf/-miss/-hwe/-r2 都应显式记录，避免框架外过滤成为黑箱。

后续 GAPIT FarmCPU/BLINK 各自独立进程和输出目录，不拼各模型最小 P，不以显著数选主模型，不把多模型一致性当独立复现。rMVP 可作为性能适配候选，不能假设不同实现效应字段天然相同。

## 11. 多重检验与结果状态

参考 harmonizer：有效 P 的个数 M，阈值 alpha/(family_size×M)。family_size=1 表示每性状控制；需全项目保守控制时填预先指定的主性状数。family_size 不替代多模型方法选择。

若部分性状失败，不能事后缩小主性状 family_size。若模型有大量未报告/数值失败测试，必须先追查原因，不能靠缩小 M 继续发表“完整主结果”。

本参考实现不计算 FDR、λGC、有效独立检验数或 permutation 阈值。将来加入需独立方法审阅。QQ/lambda 只作诊断，不自动把 lambda 调到 1；有结构/亲缘关系时普通打乱表型的 permutation 可能无效。

状态：prepared=输入已准备；software_exited_zero=外部命令退出成功；parsed=结果适配成功；partial=存在异常；complete=生产流程核对所有模块和预期产物后才能赋予。parser 成功不等于真实运行完整或科学结论有效。

## 12. 标准化关联结果

最少字段：schema_version、analysis_id、trait、model、test、genome_version、chrom、pos、variant_id、original_id、ref、alt、effect_allele、other_allele、beta、se、effect_scale、p_value、p_value_text、neg_log10_p、eaf、maf、mac、n_obs、status、filter_reason、significant。

保留 GEMMA 原始 allele1 方向，或显式转换到 ALT：反向时同时改变 beta 符号、EAF 与等位基因。等位基因/坐标不匹配直接失败，禁止静默跳过。

MAC 不能从 GEMMA 已舍入的 af 推出精确值，参考脚本填 NA。生产时加入对应性状材料的精确 allele counts。

极小 P 用 Decimal 保存，避免 1e-1000 变 0。软件已经输出 0 时标 p_underflow，不能凭空恢复精确 -log10(P)。绘图可以明确标记截断，但保存的统计表保留原始表示。

输出 gwas.tsv、significant.tsv、invalid_variants.tsv、unreported_variants.tsv、metadata.json。未报告表需要和软件日志/筛选记录对账；不能默认其中所有位点是失败或通过了检验。

参考 parser 在内存中排序，适合小/中规模验证。大规模生产应用流式处理和外部排序，按染色体分片，最终输出 BGZF/Tabix 或 Parquet。

## 13. 显著区域、LD 与候选基因

define_loci 只在选定主模型与其阈值上取 lead，LD clumping 用于去冗余而非精细定位。优先使用该性状实际分析的材料；用外部 LD 时记录参考群体匹配性。

每个区域至少登记 locus_id、trait、model、chrom、start、end、lead_variant、lead_p、boundary_method、ld_sample_sha256、参数和回退原因。无 LD 信息的固定窗口必须标识，不能称可信区间。

全群体 LD 衰减仅提供尺度参考，局部和亚群可不同。PopLDdecay/LDBlockShow 内部也有过滤默认值，需将生效值和保留标记集记录，避免其悄悄改变 SNP 集合。

annotate_loci_genes.py 执行 GFF3 gene 区间相交，输入/输出均为 1-based inclusive。flank_bp 表示在已经定义的 locus 外额外延伸，默认 0，避免重复扩大窗口。它不做变异后果注释或因果判断。

AnnoFlow 功能表、RNAFlow 表达证据另行 join；必须同组装、同基因 ID，映射多对多时保留映射证据。先做并列证据，不默认构造一个未经验证的总分。

lead genotype 表型图标注各组材料数，属于同队列描述，不是独立验证。富集背景要考虑本次可检验/可映射基因及长度/LD 偏倚。

## 14. 配置与参数映射

四件套→项目 analysisyaml→命令行最终覆盖，生产实现要做深合并和最终 schema 校验，并输出 effective_config.yaml。Genome_Version 保留框架既有命名，其余新增字段 snake_case。

parameter.qc.* 映射到 QC Shell 同名参数；parameter.structure.* 映射 pruning/PCA；traits.tsv 的协变量、变换和 n_pcs 映射 prepare；parameter.kinship.type 映射 GEMMA --kinship-type；primary_test 映射 --test；alpha/family_size 映射 harmonizer。

只有已实现参数才能进入生产 schema 的可启用分支。参考模板中列出的是目标接口，报告、local_ld 等开关仍需对应规则实现，不能直接拿模板宣称这些功能已接通。

版本、参考和资源锁在真实部署时生成：工具版本+二进制/容器摘要、Conda 明确版本、代码提交、参考/注释摘要、运行种子。不要把未验证工具版本写成已验收依赖。

## 15. Snakemake 规则参考

以下为待接入示意，不是本包已运行的 snakefile。TRAITS 上下文应来自一次解析。

~~~python
from pathlib import Path

rule prepare_trait:
    input:
        gate="01.qc/input_gate.json",
        fam="02.genotype/cohort.fam",
        samples=config["sample_csv"],
        phenotypes=config["phenotype_tsv"],
        covariates=config["covariate_tsv"],
        pcs="04.structure/pca.eigenvec"
    output:
        pheno="03.phenotype/{trait}/phenotype.gemma.txt",
        covar="03.phenotype/{trait}/covariates.gemma.txt",
        meta="03.phenotype/{trait}/metadata.json"
    params:
        script=str(Path(workflow.basedir) / "scripts" / "prepare_trait_inputs.py"),
        outdir=lambda wc: f"03.phenotype/{wc.trait}",
        n_pcs=lambda wc: TRAITS[wc.trait]["n_pcs"],
        trait=lambda wc: wc.trait
    threads: 1
    log:
        "logs/04.phenotype_qc/{trait}.log"
    shell:
        "python3 {params.script:q} "
        "--fam {input.fam:q} --samples {input.samples:q} "
        "--phenotypes {input.phenotypes:q} --trait {params.trait:q} "
        "--pcs {input.pcs:q} --n-pcs {params.n_pcs} "
        "--out-dir {params.outdir:q} > {log:q} 2>&1"
~~~

示意省略了 covariates、transform 等选项的参数构建，生产不能漏传，必须按 traits.tsv 完整构建参数数组。无协变量或 n_pcs=0 的情况不要求不存在的文件。

参考 Python 脚本拒绝覆盖非空目录；因此生产 rule 应写入本次作业的独立临时目录，成功后将明确产物移动到声明输出并清理本次临时目录。失败重试不能盲目清空整个性状或项目目录。

资源改用 rule_resource，环境按工具固定版本，所有输出包括 keep/列名/校验侧文件都要登记。以 shell command 数组构建工具命令，不允许任意自由文本参数绕过方法约束。

## 16. BioReport、交付与 FlowAgent

新增 pipeline_type=gwas 的摘要与报告适配器，报告结构包含：项目/方法、样本与位点流失、表型/结构、各性状结果、模型诊断、区域/候选证据、复现记录。

project_summary 包含模块状态、分析索引、样本/标记计数、参数、版本、警告和文件路径。报告引用实际图/表，不硬写 RNA DEG 字段。report=false 时仍可导出 JSON；GWAS 渲染器未实现前不得标 report=true 作为已支持功能。

大数据全部统计下载为文件，HTML 只加载受控规模；绘图抽样要保留全部显著位点并可审计。PDF 使用静态图。最终交付默认 copy，manifest 记录目录内文件，移走交付目录后仍可用。

FlowAgent 资产：每个实际模块 manifest、端口、最终交付目录、config_requires。端口说明补充组装、倍性、硬分型/剂量、样本顺序、性状变换和效应语义，实际兼容性由 preflight 检查，不只靠路径 pattern。

Agent 仅选择和解释已验收能力，规则库运行期只读。参数/方法修改形成可审阅计划，正式执行沿用既有提交门禁。本包不修改外部仓库或提交计算任务。

## 17. 缓存、资源与恢复

缓存身份至少包含输入/配置摘要、代码/环境、组装与注释、样本集合和顺序、实际标记集、表型变换、协变量编码、模型/检验和种子。K 还包含估计方法和 LOCO 排除集。

保持 Snakemake 默认的代码/参数/软件环境等重算依据，不直接复制只按 mtime 的示例作为可复现策略。配置、表型或协变量变化必须让对应结果失效。

每性状并行与内部 BLAS 线程共同受控。dense K 内存随 N² 增长，不能按一个固定 mem_mb 承诺所有群体。配置资源档位是起点，需由真实 benchmark 调整。

失败退出非零，按已完成合法产物恢复。允许 partial 时需要显式策略和摘要状态；不得用空 PASS、空图或空统计表掩盖严重失败。无显著结果是正常科学结果，可成功交付合法空表。

## 18. 开发工作包与验收条件

| 工作包 | 实施内容 | 验收 |
| --- | --- | --- |
| A | 输入类型扩展、四件套/schema、文件门控、身份契约 | 错误输入阻断，QC-only 可独立构建/运行 |
| B | VCF/BED 规范化、分阶段 QC、PCA、性状准备 | 顺序一致、精确身份与剔除统计 |
| C | K、GEMMA、统一结果 | 对照独立软件命令，N/效应/SE/P 一致 |
| D | 区域、LD、基因与图 | 方向、边界、组装、空结果路径正确 |
| E | 摘要、报告、复制交付 | 模块开关、无自依赖、离线交付完整 |
| F | FlowAgent、跨流程证据、LOCO/补充模型 | 单项真实验收、方法审阅与回归 |

上线前必须完成官方示例+真实植物性状+结构/亲缘相关的模拟校准；不能用包内 synthetic assoc 表代替实际运行证据。

验收覆盖：输入损坏、样本重名/顺序打乱、变换单位、协变量共线、不同性状缺失、等位基因翻转、未报告测试、零显著、部分失败、参数变更失效、报告和交付独立开关、目录产物清单和非本机路径。

## 19. 本包验证结果与限制

已通过 8 个标准库合成测试：行顺序/缺失掩码不变、重名拒绝、常量协变量拒绝、ALT 方向翻转与 1e-1000 保存、结果指纹失效、数值失败标 partial、闭区间边界/空区域表、已有输出保护。

两份 Shell 通过 bash -n；三份 CLI 支持 --help。当前没有 PLINK 2/GEMMA 二进制，因此没有执行它们，未做真实数值/性能、Snakemake、容器、集群和 BioReport 验证。

参考 Shell 不是完整审计执行器：尚未实现输入 gate、每阶段结构化剔除原因、完整命令台账、严格 K 半正定性检查、共享缓存、trait MAC、LOCO 或任意大规模流式处理。生产移植须按工作包补齐。

## 20. 官方来源

核对日期：2026-10-23。固定版本后，以该二进制 --help 和其对应源码为最终参数依据。

- [FlowFrame](https://github.com/CygnusX-Labs/FlowFrame)
- [框架规范](https://github.com/CygnusX-Labs/FlowFrame/blob/main/flow-framework/references/flow_framework_spec.md)
- [框架技能](https://github.com/CygnusX-Labs/FlowFrame/blob/main/flow-framework/SKILL.md)
- [FlowAgent 入库](https://github.com/CygnusX-Labs/FlowFrame/blob/main/docs/planning/Flow-Workflow-Onboarding-Guide.md)
- [RNAFlow](https://github.com/xsx123123/RNAFlow)
- [BCFtools 手册](https://samtools.github.io/bcftools/bcftools.html)
- [PLINK 输入](https://www.cog-genomics.org/plink/2.0/input)
- [PLINK 过滤](https://www.cog-genomics.org/plink/2.0/filter)
- [PLINK LD](https://www.cog-genomics.org/plink/2.0/ld)
- [PLINK PCA](https://www.cog-genomics.org/plink/2.0/strat)
- [GEMMA 手册](https://github.com/genetics-statistics/GEMMA/blob/master/doc/manual.tex)
- [GEMMA 参数/表型读取实现](https://github.com/genetics-statistics/GEMMA/blob/master/src/param.cpp)
- [GEMMA 选项实现](https://github.com/genetics-statistics/GEMMA/blob/master/src/gemma.cpp)
- [GAPIT](https://github.com/jiabowang/GAPIT)
- [rMVP](https://github.com/xiaolei-lab/rMVP)
- [PopLDdecay](https://github.com/BGI-shenzhen/PopLDdecay)
- [LDBlockShow](https://github.com/BGI-shenzhen/LDBlockShow)
- [Snakemake CLI](https://snakemake.readthedocs.io/en/stable/executing/cli.html)
