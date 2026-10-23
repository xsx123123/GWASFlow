# AtPolyDB → GWASFlow 测试输入

核验日期：2026-10-23。物种：Arabidopsis thaliana；参考坐标：**TAIR9**。

这个包包含原始数据下载脚本、流式转换脚本，以及已经整理好的 FT10 快速测试输入。
数据转换使用 Python 3 标准库。转成 BED/BIM/FAM 时使用 **PLINK 1.9**，命令名默认为 `plink`。
本包不包含 GWASFlow 调度器，也不会自动运行关联分析。

## 1. 直接用已整理的快速测试数据

解压这个包，进入 `AtPolyDB_GWASFlow_Test`，运行：

```bash
bash make_plink.sh data/smoke
```

若程序名为 `plink1.9`，可以运行：

```bash
PLINK_BIN=plink1.9 bash make_plink.sh data/smoke
```

转换命令会读取配套表型，并输出：

| 文件 | 用途 |
| --- | --- |
| `data/smoke/genotype.bed` | PLINK 二进制基因型 |
| `data/smoke/genotype.bim` | SNP、染色体、TAIR9 位置和观测等位基因 |
| `data/smoke/genotype.fam` | 样本顺序；FID=0，IID=ecotype ID；第6列为 FT10 |
| `data/smoke/samples.tsv` | GWASFlow 样本映射：sample、fid、iid、include |
| `data/smoke/phenotypes.tsv` | GWASFlow 表型：sample、FT10 |
| `data/smoke/phenotype.plink.tsv` | PLINK 表型：FID、IID、FT10 |
| `data/smoke/preparation.json` | 原始文件指纹、选择规则、数量和产物指纹 |

快速测试输入来自 TAIR9 原矩阵的每第20个 SNP，然后检查所选样本中的多态性与双等位性。
本次实际得到 **194 个样本、10,703 个 SNP、5 条染色体**。该稀疏子集用于检查数据流、运行和绘图；完整关联与局部 LD 测试建议用完整 SNP 集。

提供给 GWASFlow 的三个主要输入：基因型前缀 `data/smoke/genotype`、样本表 `data/smoke/samples.tsv`、表型表 `data/smoke/phenotypes.tsv`。选择连续性状 `FT10`。

## 2. 在服务器重新下载原始数据

```bash
bash download_atpolydb.sh data/raw
```

脚本下载下面两份文件，并校验本次实测 SHA256。直链已通过 HTTP HEAD 检查并实际下载：

| 文件 | 实测大小（字节） | 直链 |
| --- | --- | --- |
| `call_method_75.tar.gz` | 92,088,117 | https://raw.githubusercontent.com/Gregor-Mendel-Institute/atpolydb/master/250k_snp_data/call_method_75.tar.gz |
| `phenotype_published_raw.tsv` | 106,607 | https://raw.githubusercontent.com/Gregor-Mendel-Institute/atpolydb/master/miscellaneous_data/phenotype_published_raw.tsv |

原始地址使用上游 master，内容可能更新；脚本通过固定 SHA256 锁定本次测试快照。若日后不匹配，应先核查上游变化，勿直接跳过校验。

查看压缩包目录：

```bash
tar -tzf data/raw/call_method_75.tar.gz
```

本次实际包含：

```text
call_method_75/
call_method_75/call_method_75_info.tsv
call_method_75/call_method_75_TAIR8.csv
call_method_75/call_method_75_TAIR9.csv
```

**不必先解压两个 CSV**。转换脚本流式读取压缩包，只选 TAIR9，可避免同时展开约1.13 GB的两份矩阵。完整 TPED 另需磁盘空间。

## 3. 生成完整 FT10 测试输入

```bash
python3 prepare_atpolydb.py \
  --archive data/raw/call_method_75.tar.gz \
  --phenotypes data/raw/phenotype_published_raw.tsv \
  --trait 5_FT10 \
  --trait-name FT10 \
  --snp-step 1 \
  --out-dir data/full

bash make_plink.sh data/full
```

`data/full` 必须为新目录或空目录。转换中断后，先保留目录排查原因，再换新目录重跑。

从原始数据重新生成本包快速子集时，将 `--snp-step` 改为 `20`，并选择新的输出目录。

| 参数 | 含义 |
| --- | --- |
| `--archive` | 原始 `call_method_75.tar.gz` |
| `--phenotypes` | 原始表型 TSV |
| `--trait` | 源表型的精确列名，默认 `5_FT10` |
| `--trait-name` | 输出中的性状名，默认 `FT10`；换性状时请同时修改 |
| `--snp-step` | 每 N 个原始 SNP 选1个；1为完整集，20为快速测试 |
| `--min-samples` | 最少匹配且表型非缺失的样本数，默认20 |
| `--out-dir` | 新目录或空目录 |

## 4. 数据匹配规则

- TAIR9 矩阵有1,307个基因型样本、214,051个 SNP；染色体1–5分别为52,001、28,462、43,415、36,969、53,204个 SNP。
- 原始矩阵第一行是 array ID，第二行 `Chromosome,Positions,...` 后面才是 **ecotype ID**。
- 原始表型有199个 ecotype、107个性状；199个 ecotype 均能与基因型匹配。
- `5_FT10` 中194个样本非缺失，5个为 NA；脚本仅选择194个有效样本，并保留基因型列顺序。
- 表型原始文件含 Latin-1 字符，不能默认全部为 UTF-8。脚本先尝试 UTF-8，失败后按 Latin-1 解码；输出为 UTF-8。
- 所选基因型中的 A/C/G/T 按自交系纯合二倍体编码，例如 A→A/A；遇到未支持的编码直接报错。
- 所选样本中单态或超过两个等位基因的位点会剔除并计数；重复样本、重复坐标、非法位置会报错。
- PLINK BIM 的 A1/A2 表示观测等位基因，**没有据此断言 REF/ALT**。本包不生成带参考等位基因声明的 VCF。

## 5. 接入 GWASFlow

按项目现有配置接口映射以下值；这是输入语义，不是已经存在的统一 CLI：

| 配置含义 | 值 |
| --- | --- |
| 基因型类型 | PLINK BED/BIM/FAM |
| 基因型前缀 | `data/full/genotype`（或快速测试的 `data/smoke/genotype`） |
| 样本映射 | 同目录 `samples.tsv` |
| 表型 | 同目录 `phenotypes.tsv` |
| 性状列 / 类型 | `FT10` / quantitative |
| 物种与染色体 | 拟南芥、自交二倍体、常染色体1–5 |
| 参考坐标 | TAIR9 |
| HWE 过滤 | 自交植物测试中默认关闭 |

前一版参考包中的 `prepare_trait_inputs.py` 可以直接读取这里的样本表和表型表。`--fam` 应指向**最终QC后的FAM**，从而按最终顺序生成 GEMMA 表型和协变量。

先验收 QC、PCA、亲缘矩阵、LMM、结果解析、Manhattan/QQ 图；基因注释模块应使用配套 TAIR9 注释，不能直接混用 TAIR10 或生菜参考。抽稀子集不用于精细定位验收。

GEMMA 独立基线命令（需另装 GEMMA；该基线不加入PC协变量，与含PC的正式模型不应要求P值完全一致）：

```bash
gemma -bfile data/smoke/genotype -gk 1 -o atpolydb_FT10_K
gemma -bfile data/smoke/genotype \
  -k output/atpolydb_FT10_K.cXX.txt \
  -lmm 4 -o atpolydb_FT10_LMM
```

如需运行QC后的结果，应在以上两步共同使用同一个最终基因型文件集，并保持样本顺序一致。比较结果时固定软件版本、样本、SNP、MAF过滤、协变量、亲缘矩阵和检验类型；不能把“出现显著峰”作为唯一验收标准。

## 6. 本次验证范围

已完成：两个原始文件真实下载、SHA256、压缩包目录与表型解析、194样本×10,703 SNP真实子集生成，以及4项转换契约测试。

未完成：PLINK 1.9 二进制转换和 GEMMA 实际关联分析；本环境没有这两个程序。转换命令采用其官方接口，需在您的分析环境运行后验收。`preparation.json` 的完成状态只表示 TPED 生成完成。

运行转换契约测试：

```bash
python3 -m unittest discover -s . -p 'test_*.py' -v
```

## 来源

- 原始仓库：https://github.com/Gregor-Mendel-Institute/atpolydb
- 数据规模与参考版本：https://easygwas.biochem.mpg.de/data/public/dataset/view/1/
- 同时给出基因型与表型文件入口的研究：https://www.nature.com/articles/s41437-024-00688-z
- 原始107性状研究：Atwell et al. 2010，https://doi.org/10.1038/nature08800
- PLINK输入接口：https://www.cog-genomics.org/plink/1.9/input
- GEMMA文档：https://github.com/genetics-statistics/GEMMA
