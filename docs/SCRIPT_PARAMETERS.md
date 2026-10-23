# GWASFlow 脚本参数与软件命令手册

版本：v0.1.0-reference；日期：2026-10-23。

## 1. 阅读约定

“已提供”表示本包存在该脚本；“拟开发接口”表示开发契约，包内没有该可执行文件。不要复制拟开发命令后误认为已实现。

绝对路径示例均为占位；运行前换成实际路径。所有阈值是示例起点。命令示例不能替代参考版本确认、输入门控和真实方法审阅。

| 脚本 | 状态 | 执行环境 |
| --- | --- | --- |
| scripts/prepare_trait_inputs.py | 已实现、合成测试通过 | Python >=3.10 标准库 |
| scripts/harmonize_gemma.py | 已实现、合成测试通过 | Python >=3.10 标准库 |
| scripts/annotate_loci_genes.py | 已实现、合成测试通过 | Python >=3.10 标准库 |
| commands/01_genotype_qc.sh | 已提供、语法通过，未运行软件 | Bash、PLINK 2、Python |
| commands/02_run_gemma_reference.sh | 已提供、语法通过，未运行软件 | Bash、GEMMA、Python |
| 本文第 9 节其他名字 | 拟开发接口 | 实施时确定 |

Python 脚本错误码：0=成功；2=格式/配置/输入错误；harmonizer 另有 3=已导出结果但含数值异常。--allow-partial 可以让后者返回 0，但 metadata 仍为 partial。Shell 还会透传外部软件非零退出码，生产应统一归类。

三份 Python 输出默认保护已有文件；性状准备与标准化结果目录必须是新目录或空目录。重试应使用本次作业的新临时目录，禁止无范围地删除项目输出。

## 2. prepare_trait_inputs.py

作用：以最终 QC 的 FAM 顺序为唯一权威，对齐材料级数量性状、协变量与 PCs，建立缺失掩码和 GEMMA 输入。不读 BED，不计算 PCA，不计算 kinship。

### 2.1 全部参数

| 参数 | 类型/默认 | 必填 | 解释与校验 |
| --- | --- | --- | --- |
| --fam | 文件 | 是 | 最终 QC .fam，六列；FID/IID 二元组唯一 |
| --samples | CSV/TSV | 是 | sample,fid,iid；可选 include |
| --phenotypes | TSV | 是 | sample 和表型列；sample 唯一 |
| --trait | 字符串 | 是 | 性状内部 ID，写入 metadata |
| --phenotype-column | 字符串；同 trait | 否 | 实际表型列名 |
| --covariates | TSV | 条件 | 选择外部协变量时必须提供 |
| --numeric-covariates | 逗号列表；空 | 否 | 连续数值列，不自动标准化 |
| --categorical-covariates | 逗号列表；空 | 否 | 分类列，完整材料内做 K-1 dummy 编码 |
| --pcs | eigenvec 文件 | 条件 | n-pcs>0 时提供；接受 #FID/IID 或 #IID 表头 |
| --n-pcs | 整数；0 | 否 | 选择 PC1..PCn；不能为负 |
| --transform | none/log/log1p；none | 否 | log 要求 y>0；log1p 要求 y>-1 |
| --unit | 字符串；unspecified | 否 | 原始性状单位；效应单位还应反映变换 |
| --missing-values | 逗号列表；NA,NaN,. | 否 | 大小写不敏感；空值始终视为缺失 |
| --min-samples | 整数；20 | 否 | 最少完整材料；必须≥3，不代表足够统计功效 |
| --rank-tolerance | 数值；1e-10 | 否 | 归一化列的数值秩判定容差，0<值<1 |
| --out-dir | 新/空目录 | 是 | 输出根目录 |

数值、分类、PC 名不得重复。表型零方差、分类协变量只有一个有效水平、设计矩阵共线，或完整 N≤协变量列数+1 时阻断。含 FID 的 PCA 精确匹配二元组；没有 FID 时要求 FAM 中 IID 唯一。

未知 sample、重复 ID 和缺少 FAM 身份映射报错。表型/协变量可以包含已经在基因型 QC 中删除、但仍登记在 samples 表中的材料。

### 2.2 调用示例

~~~bash
python3 scripts/prepare_trait_inputs.py \
  --fam /project/01.workflow/qc_reference/cohort.fam \
  --samples /project/01.workflow/samples.csv \
  --phenotypes /project/01.workflow/phenotypes.tsv \
  --trait leaf_length \
  --phenotype-column leaf_length \
  --covariates /project/01.workflow/covariates.tsv \
  --categorical-covariates batch \
  --pcs /project/01.workflow/qc_reference/pca.eigenvec \
  --n-pcs 3 \
  --transform none \
  --unit cm \
  --min-samples 20 \
  --out-dir /project/01.workflow/03.phenotype/leaf_length
~~~

如果没有外部协变量，删除 --covariates 和分类/数值协变量参数；不使用 PCs 时设 --n-pcs 0 并省略 --pcs。

### 2.3 输出

| 文件 | 内容 |
| --- | --- |
| phenotype.gemma.txt | 无表头、无 ID、完整 FAM N 行；无效材料为 NA |
| covariates.gemma.txt | 无表头、无 ID、相同 N 行；第一列截距 |
| kinship_all.pheno.txt | 完整 N 行的 1，用于全队列 K |
| keep.tsv | 有效 FID/IID，无表头；供按材料提取或算 LD |
| sample_inclusion.tsv | 全部 FAM 材料和逐项排除原因 |
| sample_order.tsv | 完整 FAM 顺序的 fid/iid，有表头 |
| design_columns.tsv | 协变量矩阵列名和顺序 |
| metadata.json | 完整/有效 N、样本顺序摘要、分类编码、源文件和准备产物 SHA-256 |

这是“全矩阵+缺失掩码”路线。不要只把 phenotype 文件删成完整材料的短向量，却仍传原 N 的 BED/K，否则会破坏顺序。

## 3. 01_genotype_qc.sh

作用：从已规范化的二倍体 SNP BED 数据开始，分阶段 QC，然后生成 PCA。它不做 VCF REF/组装确认，也不建立完整 Flow 数据门控。

| 参数 | 类型/默认 | 必填 | 对应含义 |
| --- | --- | --- | --- |
| --bfile | 前缀 | 是 | 不含 .bed/.bim/.fam |
| --out-dir | 不存在的目录 | 是 | 阶段文件和日志输出目录 |
| --chr-count | 整数 | 是 | 二倍体常染色体对数，1..95；生菜示例为 9，仍需核对参考 |
| --pre-geno | 比例；0.20 | 否 | 初步宽松位点缺失阈值 |
| --mind | 比例；0.10 | 否 | 后续样本缺失阈值 |
| --geno | 比例；0.10 | 否 | 最终位点缺失阈值 |
| --maf | 比例；0.05 | 否 | 队列级最小 MAF，0..0.5 |
| --mac | 正整数；5 | 否 | 队列级最小 MAC |
| --threads | 正整数；4 | 否 | PLINK 线程上限 |
| --memory-mb | 正整数；8000 | 否 | PLINK --memory 工作内存值；实际进程还可能有额外开销 |
| --pca-components | 正整数；10 | 否 | 输出 PC 数，需小于有效矩阵秩 |
| --prune-window-kb | 正整数；50 | 否 | 剪枝物理窗口，kb |
| --prune-step | 正整数；5 | 否 | 每步向前的标记数 |
| --prune-r2 | 数值；0.20 | 否 | 剪枝 r² 阈值，(0,1] |

~~~bash
bash commands/01_genotype_qc.sh \
  --bfile /project/01.workflow/02.genotype/raw_engine \
  --out-dir /project/01.workflow/qc_reference \
  --chr-count 9 --pre-geno 0.20 --mind 0.10 --geno 0.10 \
  --maf 0.05 --mac 5 \
  --threads 4 --memory-mb 16000 \
  --pca-components 10 \
  --prune-window-kb 50 --prune-step 5 --prune-r2 0.20
~~~

内部固定 --chr-set N no-x no-y no-xy no-mt 和 --nonfounders，避免人类染色体默认值以及家系字段改变过滤人数。不运行 HWE 或亲缘删除。要求 BIM 为 1..N 数值染色体、唯一非点 ID、A/C/G/T 双等位 SNP。

顺序为 pre→sample_qc→cohort→final_stats→pruning→PCA。主输出 cohort.bed/bim/fam、final_stats.*、pca_prune.prune.in/.prune.out、pca.eigenvec/.eigenval、各阶段日志和 PLINK 版本。

## 4. 02_run_gemma_reference.sh

作用：验证准备产物摘要，计算全队列 K，然后运行一个数量性状 LMM。不是共享缓存或 LOCO 执行器。

| 参数 | 类型/默认 | 必填 | 解释 |
| --- | --- | --- | --- |
| --bfile | 最终 cohort 前缀 | 是 | 必须与准备 trait 时使用同一 FAM |
| --trait-dir | 目录 | 是 | prepare_trait_inputs 的完整输出 |
| --out-dir | 不存在的目录 | 是 | GEMMA 原始文件和运行侧信息 |
| --maf | 比例；0.05 | 否 | 显式传入 GEMMA；会按其分析材料重新判断 |
| --miss | 比例；0.10 | 否 | GEMMA 位点缺失阈值 |
| --kinship-type | centered/standardized；centered | 否 | 分别对应 -gk 1/2 |
| --test | wald/lrt/score/all；wald | 否 | 分别对应 -lmm 1/2/3/4 |
| --threads | 正整数；1 | 否 | 设置 BLAS/OpenMP 环境线程上限；不是 GEMMA 原生 --threads |

~~~bash
bash commands/02_run_gemma_reference.sh \
  --bfile /project/01.workflow/qc_reference/cohort \
  --trait-dir /project/01.workflow/03.phenotype/leaf_length \
  --out-dir /project/01.workflow/07.association/gemma/leaf_length \
  --maf 0.05 --miss 0.10 \
  --kinship-type centered --test wald --threads 1
~~~

以上 --bfile 的 FAM 与上一条 prepare 的 --fam 对应同一 qc_reference/cohort。迁移到生产结果目录时仍须保持这一约束。

内部参数映射：

| 原生参数 | 生效方式 |
| --- | --- |
| -bfile | 传入 BED 前缀 |
| -p / -n 1 | 建 K 用全 1 文件；关联用该 trait 的外部单列向量 |
| -c | 完整 N 行协变量，已含截距 |
| -gk 1/2 | centered/standardized K |
| -k | 本次全队列矩阵 |
| -lmm 1/2/3/4 | Wald/LRT/score/all |
| -maf / -miss | 显式来自 wrapper |
| -hwe 0 | 关闭 HWE |
| -r2 0.9999 | 关联时显式记录的软件协变量相关过滤阈值 |
| -outdir / -o | 目录与输出前缀，避免共享默认 output |

输出：cohort.cXX.txt 或 cohort.sXX.txt、cohort.log.txt、association.assoc.txt、association.log.txt、gemma.help.txt、run_metadata.json。

run_metadata 绑定 BED/BIM/FAM、trait metadata、K、assoc 和二进制摘要；status=software_exited_zero。仅对 K 做尺寸和有限值检查，未做 PSD、对称性或科学有效性验收。

## 5. harmonize_gemma.py

作用：核验运行/性状/assoc 身份、坐标和等位基因，导出统一统计及 Bonferroni 显著结果。

| 参数 | 类型/默认 | 必填 | 说明 |
| --- | --- | --- | --- |
| --assoc | 文件 | 是 | GEMMA 原始 assoc 表 |
| --variant-map | TSV | 是 | 必须含 variant_id,chrom,engine_chrom,pos,ref,alt |
| --sample-metadata | JSON | 是 | trait 准备输出的 metadata.json |
| --run-metadata | JSON | 是 | GEMMA wrapper 的 run_metadata.json，防结果错配 |
| --analysis-id | 字符串 | 是 | 一次方法/样本集合确定的分析 ID |
| --genome-version | 字符串 | 是 | 当前参考组装标识 |
| --effect-scale | 字符串 | 是 | 如 cm、log(cm)，不由脚本猜测 |
| --p-column | p_wald/p_lrt/p_score；p_wald | 否 | 必须属于运行记录的检验 |
| --target-effect-allele | reported/alt；reported | 否 | 保留 allele1，或变换到 ALT |
| --alpha | Decimal 字符串；0.05 | 否 | (0,1) |
| --family-size | 正整数；1 | 否 | 预声明主性状数，阈值 alpha/(family_size×有效 M) |
| --allow-partial | 开关；false | 否 | 数值异常仍标 partial，但允许返回 0 |
| --out-dir | 新/空目录 | 是 | 标准化结果目录 |

~~~bash
python3 scripts/harmonize_gemma.py \
  --assoc /project/01.workflow/07.association/gemma/leaf_length/association.assoc.txt \
  --run-metadata /project/01.workflow/07.association/gemma/leaf_length/run_metadata.json \
  --sample-metadata /project/01.workflow/03.phenotype/leaf_length/metadata.json \
  --variant-map /project/01.workflow/02.genotype/variant_map.tsv \
  --analysis-id leaf_length_gemma_primary_001 \
  --genome-version Lsat_Salinas_v11 \
  --effect-scale cm --p-column p_wald \
  --target-effect-allele reported \
  --alpha 0.05 --family-size 1 \
  --out-dir /project/01.workflow/08.results/gemma/leaf_length
~~~

原始 chr/ps 必须对应 variant_map 的 engine_chrom/pos，两个等位基因集合必须匹配。缺失映射、重复 ID、坐标或等位基因错误是结构错误，直接退出 2。

n_obs 存在则使用；否则从 n_analysis 减 n_mis/n_miss 推导并验证。af 视为报告 allele1 频率。ALT 反向时 beta×-1、eaf=1-eaf；SE 与 P 不变。

软件已输出 P=0 时标 p_underflow；保留原值，neg_log10_p 为 NA，避免编造精确度。仅从已舍入 af 无法得精确 MAC，本实现 mac=NA。

| 输出文件 | 内容 |
| --- | --- |
| gwas.tsv | 全部报告行，含异常状态 |
| significant.tsv | 当前校正阈值下的有效显著行 |
| invalid_variants.tsv | 数值无效/计数不合法行 |
| unreported_variants.tsv | variant_map 中未在 assoc 出现的 ID，待日志对账 |
| metadata.json | 有效/异常/未报告数量、阈值、指纹和状态 |

不会直接生成 lead/locus；显著位点表不等于独立关联区域表。

## 6. annotate_loci_genes.py

作用：将已经定义好的候选区域与 GFF3 gene 位置相交，按精确 gene_id 连接可选描述。

| 参数 | 类型/默认 | 必填 | 说明 |
| --- | --- | --- | --- |
| --loci | TSV | 是 | locus_id,chrom,start,end，1-based inclusive |
| --gff | 未压缩 GFF3 | 是 | 默认读取 gene feature，不支持直接 GTF |
| --genome-version | 字符串 | 是 | 声明共同组装；不会自动证明序列一致 |
| --gene-feature | 字符串；gene | 否 | GFF 第三列要匹配的特征类型 |
| --id-attribute | 字符串；ID | 否 | GFF3 属性中的唯一基因 ID |
| --flank-bp | 整数；0 | 否 | 在输入区域外再扩展的 bp，≥0 |
| --function-table | TSV | 否 | 唯一 gene_id 与 description |
| --out | 新文件路径 | 是 | 候选表；旁写 .meta.json |

~~~bash
python3 scripts/annotate_loci_genes.py \
  --loci /project/01.workflow/09.loci/leaf_length/loci.tsv \
  --gff /reference/GWASFlow_reference/Lsat_Salinas_v11/genes.gff3 \
  --genome-version Lsat_Salinas_v11 \
  --gene-feature gene --id-attribute ID \
  --flank-bp 0 \
  --function-table /reference/GWASFlow_reference/Lsat_Salinas_v11/functions.tsv \
  --out /project/01.workflow/11.candidates/leaf_length/candidates.tsv
~~~

不做 GO/KEGG 富集、转录本后果预测或因果推断。若功能表是一基因多行，先按规则汇总，不能直接 join 造成候选行重复。

输出字段：locus_id、chrom、locus_start/end、gene_id、gene_start/end、strand、distance_to_locus_bp、relation、description、status。distance 是基因与原 locus 的最短坐标距离；重叠为 0。无重叠保留 no_overlap 行；合法空 loci 输出表头。

gene chromosome 不匹配时阻断，提示检查命名/组装；当前适配器不提供“忽略匹配错误”开关。

## 7. VCF→规范化 BED 的软件命令示例

这些是独立命令示例，未封装成包内已验收的 VCF 工作流。输出文件写入前应完成参考与输入 gate；目录提前创建。本示例保留 PASS 和未评估“.”。

~~~bash
bcftools view \
  --min-alleles 2 --max-alleles 2 --types snps \
  --apply-filters 'PASS,.' \
  --output-type z --output biallelic.vcf.gz cohort.vcf.gz

bcftools norm \
  --fasta-ref genome.fa --check-ref e \
  --output-type z --output normalized.vcf.gz biallelic.vcf.gz

bcftools query \
  --format '%ID\t%CHROM\t%POS\t%REF\t%ALT\n' \
  normalized.vcf.gz > original_variant_rows.tsv

bcftools annotate \
  --set-id '%CHROM:%POS:%REF:%FIRST_ALT' \
  --output-type z --output identified.vcf.gz normalized.vcf.gz

bcftools annotate \
  --rename-chrs chromosome_rename.tsv \
  --output-type z --output engine.vcf.gz identified.vcf.gz

bcftools index --csi engine.vcf.gz

plink2 --vcf engine.vcf.gz \
  --double-id --vcf-half-call error \
  --chr-set 9 no-x no-y no-xy no-mt \
  --make-bed --out raw_engine
~~~

original_variant_rows.tsv 无表头，不是 harmonizer 可直接读取的 variant_map。生产 normalize_genotype 必须将原 ID、规范 ID、原染色体和内部染色体对应成参数手册第 5 节要求的有表头表。

chromosome_rename.tsv 是两列无表头：原参考名→内部数值名。FASTA REF 校验先在原染色体命名下进行，再重命名给引擎。不要用 --allow-extra-chr 0 把所有染色体压成 0。

FORMAT/DP、GQ 过滤另加有明确样本级语义的置缺失步骤。BCFtools expression 中针对不同样本的 AND/OR 行为需专门测试，不能把普通位点过滤表达式误当每个基因型过滤。

## 8. 原生软件参数速查

| 软件参数 | 单位/意义 | 典型误用 |
| --- | --- | --- |
| PLINK --mind | 样本缺失比例 | 和位点缺失混淆 |
| PLINK --geno | 位点缺失比例 | 假设一次调用会迭代到稳定集合 |
| PLINK --maf / --mac | 频率/计数 | 忽略所用材料集合和 founder 规则 |
| PLINK --nonfounders | 纳入全部材料统计/过滤 | 无意受 FAM 家系字段影响 |
| PLINK --chr-set | 染色体对数与类型 | 误用人类 22 条或将所有植物写 --rice |
| PLINK --indep-pairwise | 窗口、步长、r² | 将输出剪枝集当全部关联标记 |
| PLINK --pca N scols=fid | PCs，显式保留 FID | 假设输出永远有 FID |
| GEMMA -gk 1/2 | centered/standardized K | 方法不同却复用同一个缓存 |
| GEMMA -lmm 1/2/3/4 | Wald/LRT/score/all | 对每个位点取三者最小 P |
| GEMMA -p | 无 ID 表型向量 | 行顺序与 BED 不一致 |
| GEMMA -c | 无 ID 协变量矩阵 | 漏截距列 |
| GEMMA -r2 | 与协变量的相关过滤 | 误解为 LD clumping 阈值 |
| PopLDdecay -MaxDist | kb | 当成 bp |
| PopLDdecay -MAF/-Miss/-Het | 工具内部过滤 | 忽略其默认再次过滤 |
| LDBlockShow -Region | chr:start:end | 命名与输入 VCF/GFF 不一致 |
| LDBlockShow -SeleVar 2 | 以 r² 画区域 | 将 D' 与 r² 图混为一谈 |

后两个工具示例：

~~~bash
PopLDdecay -InVCF analysis_samples.vcf.gz \
  -OutStat ld_decay -MaxDist 500 -MAF 0.05 -Miss 0.10 -Het 1

LDBlockShow -InVCF analysis_samples.vcf.gz \
  -InGWAS regional_pvalues.tsv -InGFF genes.gff3 \
  -Region Chr1:100000:300000 \
  -SeleVar 2 -OutPut locus_001 -OutPng
~~~

regional_pvalues.tsv 按工具所需三列 chrom,pos,P 准备。LD 输入须使用该性状有效材料；使用外部参考时另行声明。这里的数字仅展示参数单位，不是生菜区域推荐边界。

## 9. 待开发脚本的接口契约

以下文件未提供，可作为开发任务直接拆分。CLI 名与字段是本次设计建议，不是第三方工具已存在的命令。

### 9.1 verify_inputs.py

| 参数 | 默认/约束 | 用途 |
| --- | --- | --- |
| --input-manifest | 必填 JSON | 逐文件角色、路径、组件和提供方摘要 |
| --samples / --phenotypes / --traits | 必填 | 表格契约 |
| --covariates | 可选 | 只在性状需要时必填 |
| --reference-manifest | 必填 | 组装、参考/注释摘要 |
| --expected-checksums | 可选 | 提供方校验和 |
| --checksum-mode | auto/verified_required/baseline | 无提供方基线不能冒称 verified |
| --threads | 正整数 | 校验并行度 |
| --out-dir | 新目录 | input_gate.json、md5_check.tsv、reference_check.json |

产物必须绑定生效配置和输入内容；任何不一致非零退出，禁止人工空 PASS。

### 9.2 normalize_genotype.py

| 参数 | 默认/约束 | 用途 |
| --- | --- | --- |
| --genotype-format | vcf/bcf/bed | 条件约束组件 |
| --input / --bfile | 二选一 | 基因型 |
| --gate | 必填 | 输入门控 |
| --genome-fasta / --chromosome-map | VCF 必填 | 参考与引擎染色体映射 |
| --ploidy | 2 | 其他倍性首版拒绝 |
| --representation | hardcall | 剂量输入拒绝静默转 BED |
| --multiallelic-policy | exclude | 后续 split 另设方法 |
| --filter-policy | pass_only/pass_or_unfiltered | 区分未评估 FILTER |
| --dp-min / --gq-min | null | 启用需字段存在 |
| --out-dir | 新目录 | 规范化基因型、variant_map、转换审计 |

### 9.3 compute_kinship.py / validate_kinship.py

| 参数 | 默认/约束 | 用途 |
| --- | --- | --- |
| --bfile / --sample-order | 必填 | K 与材料身份绑定 |
| --mode | genome_wide/loco | LOCO 实现前 schema 拒绝启用 |
| --exclude-chrom | loco 必填 | 受检染色体 |
| --marker-list | 可选，显式记录 | 建 K 标记集 |
| --type | centered/standardized | 估计方法 |
| --matrix / --matrix-ids | 验证器必填 | 检查行列身份 |
| --symmetry-tolerance | 1e-8，示例 | 检查前按矩阵尺度定义容差 |
| --eigenvalue-tolerance | 需方法声明 | 只容忍数值误差，不静默强制修复 |
| --out / --out-dir | 必填 | K、验证报告和侧信息 |

### 9.4 run_gapit.R

| 参数 | 默认/约束 | 用途 |
| --- | --- | --- |
| --genotype-matrix / --marker-map | 必填 | 已经固定方向与样本顺序的数值矩阵 |
| --phenotype / --covariates | 必填/可选 | 同一有效材料 |
| --sample-order | 必填 | 阻止隐式重排造成错配 |
| --trait | 必填 | 性状 |
| --model | FarmCPU 或 BLINK | 每次一个独立进程 |
| --seed | 20261023，示例 | 固定适用随机过程 |
| --threads | 正整数 | 并行上限 |
| --out-dir | 新目录 | 原始结果、运行日志、方法与环境 |

禁止在运行时 source 未固定的网络脚本；固定 GAPIT 提交与依赖，效应/P 字段做独立适配。

### 9.5 define_loci.py

| 参数 | 默认/约束 | 用途 |
| --- | --- | --- |
| --sumstats / --association-metadata | 必填 | 同一次主分析 |
| --bfile / --keep | 必填 | 该性状有效材料与标记 |
| --mode | ld_clump | 区域去冗余 |
| --p-threshold | 从 metadata 读取 | 不再另抄一份不同阈值 |
| --clump-r2 | 0.20，示例 | clump r² |
| --clump-kb | 500，示例 | 搜索半径，kb |
| --region-boundary | ld/fixed_window | 明确方法，clump 不等于精细定位 |
| --fallback-flank-bp | 显式配置 | LD 信息不足时标记回退 |
| --out-dir | 新目录 | lead.tsv、loci.tsv、locus_manifest.json |

### 9.6 run_regional_ld.py / run_ld_decay.py

| 参数 | 默认/约束 | 用途 |
| --- | --- | --- |
| --vcf / --bfile | 按后端选择 | 规范化基因型 |
| --keep | 必填 | 有效材料 |
| --loci | 区域分析必填 | 固定区域索引 |
| --max-distance-kb | 衰减分析显式配置 | 不是 bp |
| --subpopulation | 可选 | 亚群材料清单 |
| --maf / --missing-max / --het-max | 明确记录 | 同时记录 LD 工具自身过滤 |
| --metric | r2/dprime | 图和数据一致 |
| --max-regional-plots | 20，示例 | 只限制出图，不删除其余完整区域结果 |
| --out-dir | 新目录 | LD、图和逐区域状态 |

### 9.7 plot_gwas.R / plot_diagnostics.R

| 参数 | 默认/约束 | 用途 |
| --- | --- | --- |
| --sumstats / --metadata | 必填 | 统一统计与阈值 |
| --chromosome-table | 必填 | 染色体顺序和长度 |
| --trait / --model | 必填 | 图身份 |
| --pcs / --kinship / --phenotype | 诊断图按需 | 真实产物 |
| --sample-inclusion | 必填 | 图的有效材料 |
| --point-budget | 示例 200000 | 仅控制交互绘图数据量 |
| --preserve-significant | true | 显著点不得被抽样丢失 |
| --p-zero-display-cap | 显式标注 | 不修改原始统计值 |
| --formats | svg,pdf,png | 静态交付 |
| --out-dir | 新目录 | 图、绘图数据、plot_manifest.json |

### 9.8 export_project_summary.py

| 参数 | 默认/约束 | 用途 |
| --- | --- | --- |
| --effective-config | 必填 | 生效配置快照 |
| --analysis-index | 必填 | 所有主性状/模型与状态 |
| --artifact-index | 必填 | 路径、角色、摘要、来源 |
| --qc-summary / --provenance | 必填 | QC、版本与运行证据 |
| --allow-partial | false | 不把 partial 写成 complete |
| --only-qc | false | 不引用没运行的分析 |
| --schema-version | gwasflow.project.v1（提案） | 报告契约 |
| --out | 必填 | project_summary.json |

### 9.9 render_gwas_report / deliver_results

| 参数 | 默认/约束 | 用途 |
| --- | --- | --- |
| --project-summary | 渲染必填 | 统一 GWAS 摘要 |
| --pipeline-type | gwas | 不使用 RNA DEG schema |
| --format / --language | html,pdf / zh,en | 只开放已验收分支 |
| --engine / --image | 固定容器 | 按 Docker/Apptainer 环境选择 |
| --artifact-index | 交付必填 | 不用递归 glob 猜交付文件 |
| --mode | copy | 默认交付可独立移动 |
| --source-dir / --target-dir | 交付必填 | 目录范围明确 |
| --out-dir | 渲染必填 | HTML/PDF 及报告 manifest |

## 10. 配置与脚本参数对照

| 配置/表格 | 脚本参数 |
| --- | --- |
| sample_csv / phenotype_tsv / covariate_tsv | prepare --samples / --phenotypes / --covariates |
| traits.tsv: trait/phenotype_column | prepare --trait / --phenotype-column |
| traits.tsv: transform/unit | prepare --transform / --unit；harmonizer --effect-scale 需含变换 |
| traits.tsv: numeric/categorical_covariates/n_pcs | prepare 对应参数 |
| parameter.qc.pre_geno/mind/geno/maf/mac | QC Shell 对应参数 |
| parameter.structure.pca_components/prune_* | QC Shell --pca-components / --prune-* |
| parameter.kinship.type | GEMMA Shell --kinship-type |
| parameter.association.primary_test | GEMMA Shell --test |
| parameter.multiple_testing.alpha/family_size | harmonizer --alpha / --family-size |
| parameter.candidate.flank_bp | annotate --flank-bp |
| Genome_Version | harmonizer/annotate --genome-version |

--prune-step 使用标记数，--prune-window-kb、clump-kb、MaxDist 使用 kb；flank-bp 使用 bp。Schema、文档和命令构建器必须共同固定这些单位。

## 11. 最短使用顺序

1. 先阅读实现文档，确认参考、材料身份和性状设计。
2. 规范化 VCF 并生成 variant_map，或准备已核验的 BED+variant_map。
3. 运行 QC/PCA 封装，固定最终 cohort。
4. 以同一 FAM 运行 prepare_trait_inputs。
5. 运行 GEMMA 参考封装。
6. 运行 harmonizer，对账日志/未报告位点。
7. 实现 define_loci/LD 后，使用候选基因脚本。
8. 按实施工作包接入 Snakemake、FlowAgent 和 BioReport。

只有三条示例材料的 examples 表仅供格式参考，达不到 min-samples=20，且没有配套真实 BED/PCA。合成测试的完整自洽输入在 tests 中临时生成。

## 12. 官方参数来源

- [BCFtools](https://samtools.github.io/bcftools/bcftools.html)
- [PLINK 输入/染色体](https://www.cog-genomics.org/plink/2.0/input)
- [PLINK 过滤](https://www.cog-genomics.org/plink/2.0/filter)
- [PLINK LD](https://www.cog-genomics.org/plink/2.0/ld)
- [PLINK PCA](https://www.cog-genomics.org/plink/2.0/strat)
- [GEMMA 手册](https://github.com/genetics-statistics/GEMMA/blob/master/doc/manual.tex)
- [GEMMA 外部表型优先读取](https://github.com/genetics-statistics/GEMMA/blob/master/src/param.cpp)
- [GEMMA 选项实现](https://github.com/genetics-statistics/GEMMA/blob/master/src/gemma.cpp)
- [GAPIT](https://github.com/jiabowang/GAPIT)
- [PopLDdecay](https://github.com/BGI-shenzhen/PopLDdecay)
- [LDBlockShow](https://github.com/BGI-shenzhen/LDBlockShow)
