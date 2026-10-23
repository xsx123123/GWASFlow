# GWASFlow 示例模板

## 重要声明

- 示例值**仅为格式参考**：3 条材料的示例达不到 `min-samples=20`，无法运行完整分析；
  完整自洽的合成测试输入在 `tests/make_synthetic.py` 中临时生成。
- **染色体对数（9）与参考版本（Lsat_Salinas_v11）为占位值**（M2）：部署前必须替换为
  真实生菜（或目标物种）组装的对应值，并核对 `chromosome_rename.tsv` 覆盖输入 VCF
  的全部染色体。
- QC 默认阈值（mind=0.10/geno=0.10/maf=0.05/mac=5）标注为**审阅起点，非普适植物标准**；
  实际阈值须结合真实群体审阅（自交/近交植物默认关闭 HWE、不删除高亲缘材料）。

## 模板清单

| 文件 | 契约 |
| --- | --- |
| `analysisyaml.example.yaml` | 项目级配置覆盖示例（绝对路径占位） |
| `tables/samples.csv` | sample,fid,iid 必填；可选 sample_name/group/include |
| `tables/phenotypes.tsv` | 每材料一行、每性状一列；NA 缺失示范（-9 默认合法负值） |
| `tables/traits.tsv` | 九列契约全示例；两个性状示范多性状独立分析与不同 n_pcs/协变量声明 |
| `tables/covariates.tsv` | 数值（days_to_measure）+ 分类（batch，两水平）列示范 |
| `tables/chromosome_rename.tsv` | 两列无表头：原始名→数值 engine_chrom（生菜 9 对占位） |
| `tables/chromosome_table.tsv` | 染色体顺序+长度（绘图与区域展示来源；长度为示例占位） |
| `tables/expected_md5.tsv` | 提供方校验和基线示范（file/md5 两列表头） |

## traits.tsv 九列契约

`trait` / `phenotype_column`（缺省=trait）/ `trait_type`（首版仅 continuous）/
`unit` / `transform`（none|log|log1p）/ `primary_model`（首版仅 gemma_lmm）/
`numeric_covariates`（逗号列表）/ `categorical_covariates`（逗号列表）/ `n_pcs`。

性状定义单点权威：脚本参数由这张表解析生成（spec §6.3-9），不在多处重复维护。
