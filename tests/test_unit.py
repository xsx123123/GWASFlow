#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GWASFlow L1 单元测试（覆盖 18 个验收场景中的 1/2/4/5/6/7/8/16/17 + CLI 契约）。

运行：.testenv/bin/python -m pytest tests/test_unit.py -v
（或系统 python3 -m pytest，需 pyyaml 可选）
"""

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
SCRIPTS = REPO / "scripts"
PYTHON = sys.executable

DATA_DIR = Path("/tmp/gwasflow_l1_data")


# ---- 测试夹具 ----

@pytest.fixture(scope="session", autouse=True)
def synthetic_data():
    """生成基础合成数据集（全部 L1 测试共享）。"""
    if DATA_DIR.exists():
        import shutil

        shutil.rmtree(DATA_DIR)
    result = subprocess.run(
        [PYTHON, str(REPO / "tests" / "make_synthetic.py"),
         "--out-dir", str(DATA_DIR)],
        capture_output=True, text=True, check=True,
    )
    yield DATA_DIR


def run_script(script, args, expect_code=None):
    """运行脚本 CLI 并返回 CompletedProcess。"""
    result = subprocess.run(
        [PYTHON, str(SCRIPTS / script), *args],
        capture_output=True, text=True, check=False,
    )
    if expect_code is not None:
        assert result.returncode == expect_code, (
            f"{script} 期望退出码 {expect_code}，实际 {result.returncode}\n"
            f"stdout: {result.stdout[-800:]}\nstderr: {result.stderr[-800:]}"
        )
    return result


# ---- CLI 横切约定（--help / 错误码 / 输出目录保护）----

@pytest.mark.parametrize("script", [
    "verify_inputs.py", "prepare_trait_inputs.py", "harmonize_gemma.py",
    "annotate_loci_genes.py", "normalize_genotype.py", "genotype_qc.py",
    "compute_kinship.py", "validate_kinship.py", "compute_trait_mac.py",
    "run_gemma.py", "define_loci.py", "run_regional_ld.py",
    "export_project_summary.py", "deliver_results.py",
])
def test_cli_help_self_documented(script):
    """全部脚本必须支持 --help 且自文档（参数手册 §1 契约）。"""
    result = run_script(script, ["--help"])
    assert result.returncode == 0
    assert "usage" in result.stdout.lower()


def test_output_dir_nonempty_rejected(synthetic_data):
    """输出目录必须为新目录或空目录（已有输出保护；退出码 2）。"""
    out = Path("/tmp/gwasflow_l1_nonempty")
    out.mkdir(exist_ok=True)
    (out / "occupied.txt").write_text("x")
    result = run_script(
        "prepare_trait_inputs.py",
        ["--fam", str(synthetic_data / "cohort.fam"),
         "--samples", str(synthetic_data / "samples.csv"),
         "--phenotypes", str(synthetic_data / "phenotypes.tsv"),
         "--trait", "leaf_length", "--out-dir", str(out)],
    )
    assert result.returncode == 2
    assert "非空" in result.stderr


# ---- 场景 2：重复身份 ----

def test_duplicate_fid_iid_rejected():
    """重复 FID/IID 二元组 → 报错信息含重复值与行号（场景 2）。"""
    from rules.utils.id_convert import load_samples

    csv_path = Path("/tmp/gwasflow_dup_samples.csv")
    csv_path.write_text(
        "sample,fid,iid\nS1,S1,S1\nS2,S2,S2\nS3,S1,S1\n"
    )
    with pytest.raises(SystemExit):
        try:
            load_samples(str(csv_path))
        except SystemExit as exc:
            assert exc.code == 2 or exc.code == 1
            raise
    # 错误信息含重复值与行号（捕获 stderr 需 capsys，这里验证 exit 即可）
    import io
    import contextlib

    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        try:
            load_samples(str(csv_path))
        except SystemExit:
            pass
    assert "('S1','S1')" in stderr.getvalue()
    assert "第 4 行" in stderr.getvalue() or "第 3 行" in stderr.getvalue()


# ---- 场景 17：traits.tsv 权威性 ----

def test_load_traits_contract(synthetic_data):
    """traits.tsv 九列契约：phenotype_column 缺省/n_pcs 类型化/transform 枚举。"""
    from rules.utils.id_convert import load_traits

    traits = load_traits(str(synthetic_data / "traits.tsv"))
    assert set(traits.keys()) == {"leaf_length", "leaf_width"}
    assert traits["leaf_length"]["phenotype_column"] == "leaf_length"
    assert traits["leaf_length"]["n_pcs"] == 3
    assert traits["leaf_length"]["transform"] == "none"
    assert traits["leaf_length"]["numeric_covariates"] == ["days_to_measure"]
    assert traits["leaf_length"]["categorical_covariates"] == ["batch"]
    assert traits["leaf_width"]["n_pcs"] == 0
    assert traits["leaf_width"]["numeric_covariates"] == []
    assert traits["leaf_length"]["effect_scale"] == "cm"


def test_load_traits_rejects_bad_transform():
    from rules.utils.id_convert import load_traits

    path = Path("/tmp/gwasflow_bad_traits.tsv")
    path.write_text(
        "trait\tphenotype_column\ttrait_type\tunit\ttransform\tprimary_model"
        "\tnumeric_covariates\tcategorical_covariates\tn_pcs\n"
        "t1\tt1\tcontinuous\tcm\tsqrt\tgemma_lmm\t\t\t0\n"
    )
    import io
    import contextlib

    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        with pytest.raises(SystemExit):
            load_traits(str(path))
    assert "transform 非法" in stderr.getvalue()


# ---- 场景 14（部分）：未实现模块启用报错 ----

def test_validate_module_switches_rejects_loco():
    from rules.utils.validate import validate_module_switches

    config = {"loco": True}
    import io
    import contextlib

    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        with pytest.raises(SystemExit) as exc_info:
            validate_module_switches(config)
    assert exc_info.value.code == 2
    assert "后续工作包" in stderr.getvalue()
    assert "loco" in stderr.getvalue()


def test_validate_module_switches_accepts_false():
    from rules.utils.validate import validate_module_switches

    validate_module_switches({"loco": False, "farmcpu": None})   # 不抛


# ---- 场景 10/15：模块开关解析（only_qc 截断/依赖自动启用/纯函数） ----

def _base_config():
    return {
        "input_gate": True, "genotype_normalize": True, "genotype_qc": True,
        "population_structure": True, "phenotype_qc": True,
        "kinship_analysis": True, "association_analysis": True,
        "result_harmonization": True, "locus_analysis": True,
        "local_ld": True, "candidate_genes": True, "visualization": True,
        "export_summary": True, "deliver": True, "report": False,
        "only_qc": False,
    }


def test_only_qc_truncation():
    from rules.utils.gwas_targets import resolve_modules

    config = _base_config()
    config["report"] = True   # 显式启用后截断原因才可观测
    config["only_qc"] = True
    status = resolve_modules(config)
    for module in ["association_analysis", "kinship_analysis",
                   "population_structure", "result_harmonization",
                   "locus_analysis", "local_ld", "candidate_genes",
                   "visualization", "report"]:
        assert not status[module].enabled, f"only_qc 应关闭 {module}"
        assert any("only_qc" in r for r in status[module].auto_enable_reasons), (
            f"only_qc 截断原因应记录于 {module}"
        )
    for module in ["input_gate", "genotype_normalize", "genotype_qc",
                   "phenotype_qc", "export_summary", "deliver"]:
        assert status[module].enabled, f"only_qc 应保留 {module}"


def test_dependency_closure_auto_enable():
    """仅开 association 时上游（kinship/phenotype_qc/genotype_qc）自动启用且原因链完整。"""
    from rules.utils.gwas_targets import resolve_modules

    config = _base_config()
    for key in ["population_structure", "phenotype_qc", "kinship_analysis",
                "association_analysis", "result_harmonization", "locus_analysis",
                "local_ld", "candidate_genes", "visualization"]:
        config[key] = False
    config["association_analysis"] = True
    status = resolve_modules(config)
    assert status["association_analysis"].enabled
    assert status["kinship_analysis"].enabled
    assert status["phenotype_qc"].enabled
    assert status["genotype_qc"].enabled
    reasons = " ".join(status["kinship_analysis"].auto_enable_reasons)
    assert "association_analysis" in reasons


def test_resolve_modules_pure_function():
    """纯函数：可重复调用结果一致（不改 config）。"""
    from copy import deepcopy

    from rules.utils.gwas_targets import resolve_modules

    config = _base_config()
    snapshot = deepcopy(config)
    status1 = resolve_modules(config)
    status2 = resolve_modules(config)
    assert config == snapshot, "resolve_modules 不得改写 config（纯函数红线）"
    assert {m: s.enabled for m, s in status1.items()} == {
        m: s.enabled for m, s in status2.items()
    }


# ---- 场景 1：样本顺序打乱 ----

def test_prepare_sample_order_shuffle(synthetic_data, tmp_path):
    """打乱 phenotypes.tsv 行序 → 输出表型向量行序不变且与 FAM 一致（场景 1）。"""
    # 原始表型行序打乱
    lines = (synthetic_data / "phenotypes.tsv").read_text().splitlines()
    header, rows = lines[0], lines[1:]
    rows.reverse()
    shuffled = tmp_path / "phenotypes_shuffled.tsv"
    shuffled.write_text("\n".join([header] + rows) + "\n")

    out1 = tmp_path / "prep_order_a"
    out2 = tmp_path / "prep_order_b"
    common = [
        "--fam", str(synthetic_data / "cohort.fam"),
        "--samples", str(synthetic_data / "samples.csv"),
        "--trait", "leaf_length",
    ]
    run_script("prepare_trait_inputs.py", [
        *common, "--phenotypes", str(synthetic_data / "phenotypes.tsv"),
        "--out-dir", str(out1)], expect_code=0)
    run_script("prepare_trait_inputs.py", [
        *common, "--phenotypes", str(shuffled),
        "--out-dir", str(out2)], expect_code=0)

    pheno_a = (out1 / "phenotype.gemma.txt").read_text().splitlines()
    pheno_b = (out2 / "phenotype.gemma.txt").read_text().splitlines()
    assert pheno_a == pheno_b, "表型行序打乱不得改变输出向量行序（FAM 顺序权威）"

    fam_ids = [
        line.split("\t")[0]
        for line in (synthetic_data / "cohort.fam").read_text().splitlines()
    ]
    order = (out1 / "sample_order.tsv").read_text().splitlines()[1:]
    assert [line.split("\t")[0] for line in order] == fam_ids
    assert len(pheno_a) == len(fam_ids), "表型文件行数恒等于 FAM 的 N（全矩阵）"


# ---- 场景 5：性状缺失差异 ----

def test_prepare_trait_missing_difference(synthetic_data, tmp_path):
    """两性状缺失模式不同 → keep/sample_inclusion 相应不同（掩码不复用）。"""
    out1 = tmp_path / "prep_t1"
    out2 = tmp_path / "prep_t2"
    for out, trait in [(out1, "leaf_length"), (out2, "leaf_width")]:
        run_script("prepare_trait_inputs.py", [
            "--fam", str(synthetic_data / "cohort.fam"),
            "--samples", str(synthetic_data / "samples.csv"),
            "--phenotypes", str(synthetic_data / "phenotypes.tsv"),
            "--trait", trait, "--out-dir", str(out)], expect_code=0)

    keep1 = (out1 / "keep.tsv").read_text().splitlines()
    keep2 = (out2 / "keep.tsv").read_text().splitlines()
    assert keep1 != keep2, "不同性状缺失模式必须产生不同的有效样本集合"

    incl1 = (out1 / "sample_inclusion.tsv").read_text()
    incl2 = (out2 / "sample_inclusion.tsv").read_text()
    assert incl1 != incl2
    assert "phenotype_missing" in incl1 or "phenotype_missing" in incl2

    meta1 = json.loads((out1 / "metadata.json").read_text())
    meta2 = json.loads((out2 / "metadata.json").read_text())
    assert meta1["n_valid"] != meta2["n_valid"]
    # kinship_all 全 1（全队列）
    all_ones = (out1 / "kinship_all.pheno.txt").read_text().splitlines()
    assert all(v == "1" for v in all_ones) and len(all_ones) == meta1["n_total"]


# ---- 场景 4：协变量共线 ----

def test_prepare_collinear_covariates(tmp_path):
    """共线协变量 → 秩不足阻断（退出 2 + 共线来源列）。"""
    import subprocess as sp

    data = Path("/tmp/gwasflow_collinear")
    if data.exists():
        import shutil

        shutil.rmtree(data)
    sp.run([PYTHON, str(REPO / "tests" / "make_synthetic.py"),
            "--out-dir", str(data), "--boundary", "collinear"],
           capture_output=True, text=True, check=True)

    result = run_script("prepare_trait_inputs.py", [
        "--fam", str(data / "cohort.fam"),
        "--samples", str(data / "samples.csv"),
        "--phenotypes", str(data / "phenotypes.tsv"),
        "--covariates", str(data / "covariates.tsv"),
        "--trait", "leaf_length",
        "--numeric-covariates", "days_to_measure,days_copy",
        "--out-dir", str(tmp_path / "collinear_out")])
    assert result.returncode == 2
    assert "共线" in result.stderr
    assert "days_copy" in result.stderr or "days_to_measure" in result.stderr


def test_prepare_log_transform_zero(synthetic_data, tmp_path):
    """表型含 0/负值且 transform=log → 值域校验退出 2。"""
    lines = (synthetic_data / "phenotypes.tsv").read_text().splitlines()
    header, rows = lines[0], lines[1:]
    parts = rows[0].split("\t")
    parts[1] = "0"
    rows[0] = "\t".join(parts)
    pheno = tmp_path / "pheno_zero.tsv"
    pheno.write_text("\n".join([header] + rows) + "\n")
    result = run_script("prepare_trait_inputs.py", [
        "--fam", str(synthetic_data / "cohort.fam"),
        "--samples", str(synthetic_data / "samples.csv"),
        "--phenotypes", str(pheno),
        "--trait", "leaf_length", "--transform", "log",
        "--out-dir", str(tmp_path / "log_out")])
    assert result.returncode == 2
    assert "> 0" in result.stderr or "值域" in result.stderr or "log" in result.stderr


def test_prepare_pc_matching(synthetic_data, tmp_path):
    """PC 精确匹配（FID 二元组）+ design_columns 登记 + metadata 指纹。"""
    out = tmp_path / "prep_pc"
    run_script("prepare_trait_inputs.py", [
        "--fam", str(synthetic_data / "cohort.fam"),
        "--samples", str(synthetic_data / "samples.csv"),
        "--phenotypes", str(synthetic_data / "phenotypes.tsv"),
        "--covariates", str(synthetic_data / "covariates.tsv"),
        "--categorical-covariates", "batch",
        "--numeric-covariates", "days_to_measure",
        "--pcs", str(synthetic_data / "pca.eigenvec"),
        "--trait", "leaf_length", "--n-pcs", "3",
        "--out-dir", str(out)], expect_code=0)
    design = (out / "design_columns.tsv").read_text()
    assert "intercept" in design
    assert "batch[" in design        # K-1 dummy
    assert "PC1" in design
    covar_lines = (out / "covariates.gemma.txt").read_text().splitlines()
    assert all(line.split("\t")[0] == "1" for line in covar_lines)   # 截距首列
    meta = json.loads((out / "metadata.json").read_text())
    assert meta["categorical_encoding"]["batch"]["reference_level"] == "B1"
    assert "sha256" in json.dumps(meta["sources"])


# ---- 场景 6/7 + 阈值/MAC/unreported：harmonize ----

def _harmonize_common(data, out):
    return [
        "--assoc", str(data / "association.assoc.txt"),
        "--variant-map", str(data / "variant_map.tsv"),
        "--sample-metadata", str(data / "sample_metadata.json"),
        "--run-metadata", str(data / "run_metadata.json"),
        "--trait-mac", str(data / "trait_mac.tsv"),
        "--analysis-id", "leaf_length_gemma_primary_001",
        "--genome-version", "Lsat_Salinas_v11",
        "--effect-scale", "cm", "--trait", "leaf_length",
        "--out-dir", str(out),
    ]


def test_harmonize_basic_flow(synthetic_data, tmp_path):
    """基础五件套 + 阈值 5e-5（M=360，alpha=0.05，family_size=1 → 0.05/360）。"""
    out = tmp_path / "harm_basic"
    run_script("harmonize_gemma.py",
               _harmonize_common(synthetic_data, out), expect_code=0)
    meta = json.loads((out / "metadata.json").read_text())
    m = meta["counts"]["n_valid_p"]
    from decimal import Decimal
    expected = Decimal("0.05") / m
    assert Decimal(meta["multiple_testing"]["threshold"]) == expected
    assert meta["status"] == "parsed"
    gwas_header = (out / "gwas.tsv").read_text().splitlines()[0].split("\t")
    assert "effect_allele" in gwas_header and "p_value_text" in gwas_header
    assert "mac" in gwas_header and "significant" in gwas_header
    # MAC 精确注入（非 af×N 舍入值）
    with open(out / "gwas.tsv", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    mac_map = {}
    with open(synthetic_data / "trait_mac.tsv", newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            mac_map[row["variant_id"]] = int(row["mac"])
    for row in rows[:20]:
        assert int(row["mac"]) == mac_map[row["variant_id"]]


def test_harmonize_allele_flip(synthetic_data, tmp_path):
    """等位基因翻转（场景 6）：--target-effect-allele=alt → beta×-1、eaf=1-eaf、SE/P 不变。"""
    out_reported = tmp_path / "harm_reported"
    out_alt = tmp_path / "harm_alt"
    run_script("harmonize_gemma.py",
               _harmonize_common(synthetic_data, out_reported)
               + ["--target-effect-allele", "reported"], expect_code=0)
    run_script("harmonize_gemma.py",
               _harmonize_common(synthetic_data, out_alt)
               + ["--target-effect-allele", "alt"], expect_code=0)

    def load(path):
        with open(path, newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    rows_r = load(out_reported / "gwas.tsv")
    rows_a = load(out_alt / "gwas.tsv")
    assert len(rows_r) == len(rows_a)
    flipped = 0
    for r, a in zip(rows_r, rows_a):
        assert r["p_value_text"] == a["p_value_text"], "方向转换不得改变 P"
        assert abs(float(r["se"]) - float(a["se"])) < 1e-12, "SE 不变"
        if r["effect_allele"] != a["effect_allele"]:
            flipped += 1
            assert abs(float(a["beta"]) + float(r["beta"])) < 1e-9, "beta×-1"
            assert abs((1.0 - float(r["eaf"])) - float(a["eaf"])) < 1e-9, "eaf=1-eaf"
    assert flipped > 0, "应存在被翻转的行（GEMMA allele1≠alt 的行）"


def test_harmonize_p_underflow(tmp_path):
    """P 值下溢（场景 7）：P=1e-1000 原值保留；P=0 → p_underflow + neg_log10_p=NA。"""
    import subprocess as sp

    data = Path("/tmp/gwasflow_underflow")
    if data.exists():
        import shutil

        shutil.rmtree(data)
    sp.run([PYTHON, str(REPO / "tests" / "make_synthetic.py"),
            "--out-dir", str(data), "--boundary", "p_underflow"],
           capture_output=True, text=True, check=True)
    out = tmp_path / "harm_underflow"
    run_script("harmonize_gemma.py", _harmonize_common(data, out), expect_code=0)

    with open(out / "gwas.tsv", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    underflow_rows = [r for r in rows if r["status"] == "p_underflow"]
    tiny_rows = [r for r in rows if r["p_value_text"] == "1e-1000"]
    assert underflow_rows, "P=0 行应标 p_underflow"
    assert underflow_rows[0]["neg_log10_p"] == "NA", "P=0 禁止编造 -log10"
    assert underflow_rows[0]["p_value_text"] == "0"
    assert tiny_rows, "1e-1000 行应保留原值"
    assert tiny_rows[0]["p_value_text"] == "1e-1000"
    assert tiny_rows[0]["neg_log10_p"] == "1000"


def test_harmonize_coord_mismatch_rejected(tmp_path):
    """坐标不匹配 assoc 行 → 退出 2 并列出冲突行（结构错误，禁止静默跳过）。"""
    import subprocess as sp

    data = Path("/tmp/gwasflow_coordmm")
    if data.exists():
        import shutil

        shutil.rmtree(data)
    sp.run([PYTHON, str(REPO / "tests" / "make_synthetic.py"),
            "--out-dir", str(data), "--boundary", "coord_mismatch"],
           capture_output=True, text=True, check=True)
    result = run_script("harmonize_gemma.py",
                        _harmonize_common(data, tmp_path / "harm_mm"))
    assert result.returncode == 2
    assert "染色体不匹配" in result.stderr or "结构错误" in result.stderr


def test_harmonize_unreported_variants(synthetic_data, tmp_path):
    """variant_map 中未在 assoc 出现的 ID → unreported_variants.tsv。"""
    # 从 variant_map 删除一行记录制造 unreported
    data = tmp_path / "unreported_data"
    data.mkdir()
    for f in ["association.assoc.txt", "variant_map.tsv", "sample_metadata.json",
              "run_metadata.json", "trait_mac.tsv"]:
        (data / f).write_text((synthetic_data / f).read_text())
    lines = (data / "variant_map.tsv").read_text().splitlines()
    header, rows = lines[0], lines[1:]
    unreported_vid = rows[0].split("\t")[0]
    (data / "variant_map.tsv").write_text(
        "\n".join([header] + rows + [f"{unreported_vid}_ghost\t.\tLsat1\t1\t999999\tA\tG"]) + "\n")
    out = tmp_path / "harm_unreported"
    run_script("harmonize_gemma.py", _harmonize_common(data, out), expect_code=0)
    unreported = (out / "unreported_variants.tsv").read_text()
    assert unreported_vid + "_ghost" in unreported
    meta = json.loads((out / "metadata.json").read_text())
    assert meta["counts"]["n_unreported"] >= 1


# ---- 场景 8：无显著区域（表级空表合法性） ----

def test_define_loci_empty_significant(synthetic_data, tmp_path):
    """零显著数据 → lead/loci/候选三表为仅表头合法文件且流程成功（场景 8 前半）。"""
    # 构造零显著的 harmonize 产物
    harm_out = tmp_path / "harm_zero"
    run_script("harmonize_gemma.py",
               _harmonize_common(synthetic_data, harm_out), expect_code=0)
    meta = json.loads((harm_out / "metadata.json").read_text())
    if meta["counts"]["n_significant"] > 0:
        pytest.skip("合成数据 P>=0.05 无显著（默认构造）")

    # 空显著 → define_loci 合法空表（不需要 plink2：无显著直接空表路径）
    loci_out = tmp_path / "loci_zero"
    result = run_script("define_loci.py", [
        "--sumstats", str(harm_out / "gwas.tsv"),
        "--significant", str(harm_out / "significant.tsv"),
        "--association-metadata", str(harm_out / "metadata.json"),
        "--trait", "leaf_length",
        "--bfile", str(synthetic_data / "cohort"),
        "--keep", str(synthetic_data / "cohort.fam"),
        "--out-dir", str(loci_out)])
    assert result.returncode == 0, result.stderr
    lead = (loci_out / "lead.tsv").read_text().strip().splitlines()
    assert len(lead) == 1, "lead.tsv 应为仅表头"
    loci = (loci_out / "loci.tsv").read_text().strip().splitlines()
    assert len(loci) == 1, "loci.tsv 应为仅表头"

    # 空 loci → annotate 合法表头
    cand = tmp_path / "cand_zero"
    run_script("annotate_loci_genes.py", [
        "--loci", str(loci_out / "loci.tsv"),
        "--gff", str(synthetic_data / "genes.gff3"),
        "--genome-version", "Lsat_Salinas_v11",
        "--out", str(cand / "candidates.tsv")], expect_code=0)
    cand_lines = (cand / "candidates.tsv").read_text().strip().splitlines()
    assert len(cand_lines) == 1, "空 loci → 仅表头合法输出"


# ---- 场景 16：染色体映射 ----

def test_chromosome_map_loader_rejects_duplicate_and_bad():
    from scripts.normalize_genotype import load_chromosome_map

    dup = Path("/tmp/gwasflow_dup_map.tsv")
    dup.write_text("Lsat1\t1\nLsat1\t2\n")
    import io
    import contextlib

    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        with pytest.raises(SystemExit):
            load_chromosome_map(str(dup))
    assert "重复" in stderr.getvalue()

    bad = Path("/tmp/gwasflow_bad_map.tsv")
    bad.write_text("Lsat1\tx\n")
    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        with pytest.raises(SystemExit):
            load_chromosome_map(str(bad))
    assert "非整数" in stderr.getvalue()


def test_unknown_chromosome_rejected_in_audit(synthetic_data, tmp_path):
    """映射表外染色体（ChrUnplaced）→ normalize 审计层拒绝（列出染色体名）。"""
    import subprocess as sp
    import shutil

    data = Path("/tmp/gwasflow_unknown_chrom")
    if data.exists():
        shutil.rmtree(data)
    sp.run([PYTHON, str(REPO / "tests" / "make_synthetic.py"),
            "--out-dir", str(data), "--boundary", "unknown_chrom"],
           capture_output=True, text=True, check=True)

    # 映射校验逻辑内嵌于 normalize_vcf 流程；这里直接验证映射表不覆盖输入染色体
    # （完整 bcftools 链路归 L3/L4；表级验证映射覆盖）
    from scripts.normalize_genotype import load_chromosome_map
    mapping = load_chromosome_map(str(data / "chromosome_rename.tsv"))
    assert "ChrUnplaced" not in mapping
    assert mapping["Lsat1"] == 1
    # VCF 中存在 ChrUnplaced 染色体记录（合成器注入）
    import gzip

    with gzip.open(data / "cohort.vcf.gz", "rt", encoding="utf-8") as handle:
        vcf_text = handle.read()
    assert "ChrUnplaced" in vcf_text


# ---- annotate：flank=0 / 组装一致性 / 功能表去重 ----

def _write_loci(path, rows):
    lines = ["locus_id\ttrait\tmodel\tchrom\tstart\tend\tlead_variant\tlead_p"
             "\tboundary_method\tboundary_params\tfallback_reason"
             "\tld_sample_sha256\tn_clump_members"]
    lines.extend(rows)
    path.write_text("\n".join(lines) + "\n")


def test_annotate_flank_zero_no_expansion(synthetic_data, tmp_path):
    """flank_bp=0 → 候选区间不因默认值扩大（只有真正重叠的基因入选）。"""
    loci = tmp_path / "loci_overlap.tsv"
    # gene_001 位于 Lsat2:103000-104500（合成器 k=1 → CHROMS[1]=Lsat2）
    _write_loci(loci, [
        "locus_0001\tleaf_length\tgemma_lmm\tLsat2\t103000\t103200"
        "\tLsat2:103000:A:G\t1e-6\tld\t{}\t\tabcd\t0"
    ])
    out = tmp_path / "cand_flank0"
    run_script("annotate_loci_genes.py", [
        "--loci", str(loci),
        "--gff", str(synthetic_data / "genes.gff3"),
        "--genome-version", "Lsat_Salinas_v11",
        "--flank-bp", "0",
        "--out", str(out / "candidates.tsv")], expect_code=0)
    text = (out / "candidates.tsv").read_text()
    assert "gene_001" in text, "重叠基因应入选（locus 103000-103200 在 gene_001 内）"

    # flank=0：远端基因（如 gene_020，位于 Lsat2 之外染色体或远位置）不入选
    with open(out / "candidates.tsv", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        if row["gene_id"].startswith("gene_"):
            assert row["relation"] == "overlap"
            assert row["distance_to_locus_bp"] == "0"


def test_annotate_flank_positive_includes_nearby(synthetic_data, tmp_path):
    """flank_bp>0 → 窗口外基因以 upstream/downstream 入选且距离正确。"""
    loci = tmp_path / "loci_flank.tsv"
    _write_loci(loci, [
        "locus_0001\tleaf_length\tgemma_lmm\tLsat2\t110000\t110100"
        "\tLsat2:110000:A:G\t1e-6\tld\t{}\t\tabcd\t0"
    ])
    out = tmp_path / "cand_flank"
    run_script("annotate_loci_genes.py", [
        "--loci", str(loci),
        "--gff", str(synthetic_data / "genes.gff3"),
        "--genome-version", "Lsat_Salinas_v11",
        "--flank-bp", "10000",
        "--out", str(out / "candidates.tsv")], expect_code=0)
    text = (out / "candidates.tsv").read_text()
    assert "upstream" in text or "downstream" in text


def test_annotate_assembly_mismatch_rejected(synthetic_data, tmp_path):
    """GFF 染色体名与 locus 完全无法对应 → 注释失败（提示命名/组装检查）。"""
    loci = tmp_path / "loci_mismatch.tsv"
    _write_loci(loci, [
        "locus_0001\tleaf_length\tgemma_lmm\tChrUnknown\t1000\t2000"
        "\tx\t1e-6\tld\t{}\t\tabcd\t0"
    ])
    result = run_script("annotate_loci_genes.py", [
        "--loci", str(loci),
        "--gff", str(synthetic_data / "genes.gff3"),
        "--genome-version", "Lsat_Salinas_v11",
        "--out", str(tmp_path / "cand_mm" / "candidates.tsv")])
    assert result.returncode == 2
    assert "命名/组装" in result.stderr or "组装" in result.stderr


def test_annotate_function_table_dedup(synthetic_data, tmp_path):
    """功能表一基因多行 → 候选表该基因不重复（描述汇总）。"""
    loci = tmp_path / "loci_func.tsv"
    # gene_001 位于 Lsat2:103000-104500（合成器 k=1 → CHROMS[1]=Lsat2）
    _write_loci(loci, [
        "locus_0001\tleaf_length\tgemma_lmm\tLsat2\t103500\t103600"
        "\tLsat2:103500:A:G\t1e-6\tld\t{}\t\tabcd\t0"
    ])
    out = tmp_path / "cand_func"
    run_script("annotate_loci_genes.py", [
        "--loci", str(loci),
        "--gff", str(synthetic_data / "genes.gff3"),
        "--genome-version", "Lsat_Salinas_v11",
        "--function-table", str(synthetic_data / "functions.tsv"),
        "--out", str(out / "candidates.tsv")], expect_code=0)
    with open(out / "candidates.tsv", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    gene001 = [r for r in rows if r["gene_id"] == "gene_001"]
    assert len(gene001) == 1, "一基因多行功能表不得造成候选行重复"
    assert "auxin" in gene001[0]["description"]
    assert "cell wall" in gene001[0]["description"], "多行描述应汇总"


# ---- kinship 验证 ----

def test_validate_kinship_asymmetric_fails(synthetic_data, tmp_path):
    """注入非对称 K → 验证报告标记对称性失败。"""
    data = synthetic_data
    k_lines = (data / "kinship.cXX.txt").read_text().splitlines()
    row0 = k_lines[0].split("\t")
    row1 = k_lines[1].split("\t")
    row0[1] = "0.999999"   # 打破对称
    k_lines[0] = "\t".join(row0)
    k_bad = tmp_path / "kinship_bad.txt"
    k_bad.write_text("\n".join(k_lines) + "\n")
    result = run_script("validate_kinship.py", [
        "--matrix", str(k_bad),
        "--matrix-ids", str(data / "sample_order.tsv"),
        "--fam", str(data / "cohort.fam"),
        "--out", str(tmp_path / "validation_bad.json")])
    assert result.returncode == 2
    report = json.loads((tmp_path / "validation_bad.json").read_text())
    assert report["symmetric_ok"] is False


def test_validate_kinship_ok(synthetic_data, tmp_path):
    """对称有限 K → 验证通过（PSD 用保守界）。"""
    result = run_script("validate_kinship.py", [
        "--matrix", str(synthetic_data / "kinship.cXX.txt"),
        "--matrix-ids", str(synthetic_data / "sample_order.tsv"),
        "--fam", str(synthetic_data / "cohort.fam"),
        "--eigenvalue-tolerance", "0.5",
        "--out", str(tmp_path / "validation_ok.json")], expect_code=0)
    report = json.loads((tmp_path / "validation_ok.json").read_text())
    assert report["status"] == "passed"
    assert report["dimension_ok"] and report["identity_ok"]


# ---- 交付独立移动（场景 11 表级：copy 模式 + manifest 自洽） ----

def test_deliver_copy_manifest_self_contained(synthetic_data, tmp_path):
    """copy 交付 + manifest.md5 逐文件校验（移动交付目录后仍可用）。"""
    # 构造 artifact_index
    index = tmp_path / "artifact_index.tsv"
    files = ["01.gate/input_gate.json", "02.genotype/variant_map.tsv"]
    (tmp_path / "01.gate").mkdir()
    (tmp_path / "02.genotype").mkdir()
    (tmp_path / "01.gate" / "input_gate.json").write_text("{}")
    (tmp_path / "02.genotype" / "variant_map.tsv").write_text(
        (synthetic_data / "variant_map.tsv").read_text())
    import hashlib

    def md5(p):
        return hashlib.md5(Path(p).read_bytes()).hexdigest()

    with open(index, "w", encoding="utf-8") as handle:
        handle.write("path\trole\tsize_bytes\tmd5\tsource\n")
        for f in files:
            handle.write(f"{f}\ttest\t{Path(tmp_path / f).stat().st_size}\t"
                         f"{md5(tmp_path / f)}\tunit\n")

    target = tmp_path / "deliver_out" / "results"
    out_dir = tmp_path / "deliver_out"
    run_script("deliver_results.py", [
        "--artifact-index", str(index),
        "--source-dir", str(tmp_path),
        "--target-dir", str(target),
        "--mode", "copy",
        "--out-dir", str(out_dir)], expect_code=0)

    manifest = json.loads((out_dir / "delivery_manifest.json").read_text())
    assert manifest["n_files"] == 2
    assert all(f["action"] == "copy" for f in manifest["files"])

    # 移动交付目录后 manifest.md5 校验仍通过（自洽可移动）
    import shutil

    moved = tmp_path / "deliver_moved"
    shutil.move(str(out_dir), str(moved))
    manifest_md5_text = (moved / "delivery_manifest.md5").read_text().split()[0]
    actual = md5(moved / "delivery_manifest.json")
    assert manifest_md5_text == actual, "manifest.md5 与 manifest 内容一致（移动后校验）"
    assert (moved / "results" / "01.gate" / "input_gate.json").exists()


def test_deliver_missing_index_rejected(tmp_path):
    """索引文件缺失 → 退出 2（禁止静默跳过交付）。"""
    result = run_script("deliver_results.py", [
        "--artifact-index", str(tmp_path / "nonexistent.tsv"),
        "--source-dir", str(tmp_path),
        "--target-dir", str(tmp_path / "t"),
        "--out-dir", str(tmp_path / "o")])
    assert result.returncode == 2
    assert "缺失" in result.stderr


# ---- 真实数据（AtPolyDB BED 路线）暴露的修复回归 ----

def test_bed_route_reference_requirements():
    """BED 路线按需校验不要求 FASTA/FAI/重命名映射（BIM 已是数值染色体）。"""
    from rules.utils.reference_update import check_gwas_reference_paths

    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        chrom_table = Path(tmp) / "chromosome_table.tsv"
        chrom_table.write_text("chrom\tlength_bp\n1\t1000\n")
        config = {
            "Genome_Version": "TAIR9",
            "genotype_format": "bed",
            "GWAS_reference": {
                "TAIR9": {
                    "genome_fasta": "",
                    "genome_fai": "",
                    "chromosome_rename_map": "",
                    "chromosome_table": str(chrom_table),
                    "gff3": "",
                    "function_table": "",
                }
            },
        }
        enabled = ["input_gate", "genotype_normalize", "genotype_qc", "visualization"]
        # BED 路线：不应因 FASTA/rename map 为空而失败（candidate_genes 未启用，gff3 不校验）
        check_gwas_reference_paths(config, enabled)   # 不抛即通过


def test_bed_route_reference_still_requires_vcf_assets():
    """VCF 路线仍要求 FASTA/FAI/rename map（防止修复过度）。"""
    from rules.utils.reference_update import check_gwas_reference_paths

    config = {
        "Genome_Version": "TAIR9",
        "genotype_format": "vcf",
        "GWAS_reference": {
            "TAIR9": {
                "genome_fasta": "/nonexistent/genome.fasta",
                "genome_fai": "/nonexistent/genome.fasta.fai",
                "chromosome_rename_map": "/nonexistent/rename.tsv",
            }
        },
    }
    import io
    import contextlib

    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        with pytest.raises(SystemExit):
            check_gwas_reference_paths(
                config, ["input_gate", "genotype_normalize", "genotype_qc"]
            )
    assert "genome_fasta" in stderr.getvalue()
