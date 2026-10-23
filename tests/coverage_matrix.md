# GWASFlow 验收场景覆盖矩阵（18 场景 ↔ 测试用例 ↔ spec 规则）

| # | 验收场景 | L1 单元测试 | L2 DAG 测试 | L3 行为测试 | 实现落点 |
| --- | --- | --- | --- | --- | --- |
| 1 | 样本顺序打乱 | test_prepare_sample_order_shuffle | —（FAM 顺序权威为单元级） | 端到端 FAM 对齐 | prepare_trait_inputs.py §2.1.3-3 |
| 2 | 重复身份 | test_duplicate_fid_iid_rejected | — | 门控阻断（端到端含 samples 校验） | id_convert.py + verify_inputs.py |
| 3 | 输入内容变化 | — | 场景 3（重算链） | 失败恢复（QC 产物保留） | verify_inputs.py 门控绑定 |
| 4 | 协变量共线 | test_prepare_collinear_covariates | — | — | prepare_trait_inputs.py 秩检查 |
| 5 | 性状缺失差异 | test_prepare_trait_missing_difference | 两性状独立展开 | — | prepare + {trait} wildcard |
| 6 | 等位基因翻转 | test_harmonize_allele_flip | — | — | harmonize_gemma.py 方向转换 |
| 7 | P 值下溢 | test_harmonize_p_underflow | — | — | harmonize_gemma.py Decimal |
| 8 | 无显著区域 | test_define_loci_empty_significant + annotate 空表 | — | 端到端空表路径 | define_loci.py + annotate |
| 9 | 部分失败 | —（脚本级 failed 传导见元数据校验） | — | 场景 9 三断言（failed_software/partial/健康性状） | run_gemma.py --keep-going + harmonize 传导 |
| 10 | 模块开关 | test_only_qc_truncation + test_dependency_closure_auto_enable + 纯函数 | local_ld/candidate_genes 关闭断言 | — | gwas_targets.py resolve_modules |
| 11 | 交付独立移动 | test_deliver_copy_manifest_self_contained | — | 场景 11 两断言 | deliver_results.py copy 模式 |
| 12 | 参数变化后重算 | — | 场景 12（n_pcs 变化） | 输入变化重跑 | Snakemake 重算依据（无 mtime 缓存） |
| 13 | GEMMA 数值对账 | —（未验证项 M3） | — | —（mock 不做数值断言） | run_gemma.py D10 bfile 全集可比 |
| 14 | 未实现模块启用报错 | test_validate_module_switches_rejects_loco | 场景 14（Schema 层拒绝） | — | schema enum + validate_module_switches |
| 15 | only_qc 模式 | test_only_qc_truncation | 场景 15（DAG 截断五断言） | — | resolve_modules 截断 + skip-pcs 降级 |
| 16 | 染色体映射 | test_chromosome_map_loader_rejects_* + unknown_chrom 审计 | — | —（bcftools mock 透传 rename） | normalize_genotype.py 映射校验 |
| 17 | traits.tsv 权威性 | test_load_traits_contract + bad_transform | 端到端参数派生 | — | id_convert.py load_traits 单点 |
| 18 | 软件过滤对账 | —（未验证项 M3） | — | run_metadata 四过滤记录（产物断言） | run_gemma.py software_filters |

## 汇总

- L1：44 项断言通过（pytest）
- L2：21 项断言通过（snakemake dry-run）
- L3：25 项断言通过（mock 外部软件端到端）
- L4/L5：**未验证项**（M3/M4 缺真实软件与真实数据；记录进 project_summary.json
  unverified_items，不冒称已验证）
