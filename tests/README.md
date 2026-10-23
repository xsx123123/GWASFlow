# GWASFlow 分层测试体系

## 五层测试（design §2.1.3-8；FRAMEWORK_ALIGNMENT.md §6 分层验收策略）

| 层 | 入口 | 覆盖 | 外部依赖 |
| --- | --- | --- | --- |
| L1 单元 | `.testenv/bin/python -m pytest tests/test_unit.py -v` | 场景 1/2/4/5/6/7/8/16/17 + CLI 契约（--help/错误码/输出目录保护） | 无 |
| L2 DAG | `bash tests/test_dag.sh` | 场景 3/10/12/14/15 + 两性状无互依赖边、基线无循环依赖 | snakemake（.testenv） |
| L3 行为 | `bash tests/test_behavior.sh` | 场景 9/11 + 门控三态、失败恢复、effective_config 生效、全链产物齐备 | mock 工具（自动生成） |
| L4 外部软件对账 | `tests/commands/`（M3 具备时执行） | 场景 13（GEMMA 数值对账 N/beta/SE/P/标记集）、18（软件过滤对账） | 真实 PLINK2/GEMMA/BCFtools |
| L5 端到端真实数据 | —（M4 具备时执行） | 真实植物群体上线前校准 | 真实数据 + 集群 |

## 运行环境准备

```bash
# 测试虚拟环境（snakemake + pytest）
uv venv .testenv --python 3.12
uv pip install --python .testenv/bin/python pytest pyyaml "snakemake>=9.9.0"

# L1
.testenv/bin/python -m pytest tests/test_unit.py -v

# L2（Snakemake dry-run DAG 断言）
bash tests/test_dag.sh

# L3（mock plink2/gemma/bcftools/Rscript 驱动真实规则链端到端）
bash tests/test_behavior.sh
```

## 测试资产

| 文件 | 用途 |
| --- | --- |
| `make_synthetic.py` | 合成数据生成器（小群体 BED/VCF/表型/协变量/GFF3/assoc/K 等）+ 边界构造器（--boundary） |
| `make_mock_tools.py` | L3 mock 外部软件生成器（plink2/gemma/bcftools/Rscript） |
| `test_unit.py` | L1 pytest 用例（44 项） |
| `test_dag.sh` | L2 DAG 断言（21 项） |
| `test_behavior.sh` | L3 端到端行为（25 项） |
| `coverage_matrix.md` | 18 个验收场景 ↔ 测试用例 ↔ spec 规则三方映射矩阵 |

## 外部软件就绪检测与未验证项

- L4 层执行前检测：`command -v plink2 && command -v gemma && command -v bcftools`
  （版本固定后以 `--version` 输出与 `tests/commands/` 独立命令对账）
- L4/L5 不具备时**记录为未验证项**（进 project_summary.json 的 `unverified_items`，
  不冒称已验证）：
  - M3：真实软件数值对账（场景 13/18）
  - M4：真实群体端到端
  - M5：rich-loguru 插件生产接线
  - M6：BioReport 渲染器验收

## 契约

- 合成数据与 mock 工具只用于相应测试，**不能替代真实软件验收**（M3）。
- 覆盖证明见 `tests/coverage_matrix.md`（18 场景无遗漏映射）。
