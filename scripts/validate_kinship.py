#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GWASFlow kinship 矩阵验证（参数手册 §9.3）。

检查项：
1. N×N 维度（与 sample_order/FAM 一致）
2. 对称性（容差按矩阵尺度声明：tol × max(1, |K|max)）
3. 行列身份（sample_order 与 FAM 顺序逐项一致）
4. 数值有限性
5. PSD 特征值检查（最小特征值 >= -tol_eig × max(1, |eig|max)；仅容忍数值误差级别，
   禁止静默强制修复）

错误码：0=验证通过；2=验证失败（维度/顺序/有限性/对称性任一失败）。
PSD 仅容忍数值误差：超容忍即失败退出 2。
"""

import argparse
import json
import sys
from pathlib import Path

EXIT_OK = 0
EXIT_INPUT_ERROR = 2


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="validate_kinship.py",
        description="kinship 矩阵维度/对称性/身份/有限性/PSD 检查（不静默修复）。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--matrix", required=True, help="K 矩阵文件（空白分隔 N×N）")
    parser.add_argument("--matrix-ids", required=True, help="sample_order.tsv（fid/iid 有表头）")
    parser.add_argument("--fam", required=True, help="cohort.fam（顺序对账基准）")
    parser.add_argument("--symmetry-tolerance", type=float, default=1e-8,
                        help="对称性相对容差（按矩阵尺度缩放）")
    parser.add_argument("--eigenvalue-tolerance", type=float, default=1e-8,
                        help="PSD 最小特征值相对容忍（仅数值误差级别）")
    parser.add_argument("--out", required=True, help="验证报告 JSON 输出路径")
    return parser.parse_args(argv)


def fail(message):
    print(f"[validate_kinship ERROR] {message}", file=sys.stderr)
    sys.exit(EXIT_INPUT_ERROR)


def jacobi_min_eigenvalue(matrix, max_sweeps=100, tol=1e-10):
    """对称矩阵最小特征值（Jacobi 旋转，标准库实现，用于无 numpy 环境）。"""
    import copy
    import math

    n = len(matrix)
    a = [row[:] for row in matrix]
    for _ in range(max_sweeps):
        off = math.sqrt(sum(a[i][j] ** 2 for i in range(n) for j in range(i + 1, n)))
        if off < tol:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                if abs(a[p][q]) < tol * 1e-6:
                    continue
                theta = (a[q][q] - a[p][p]) / (2 * a[p][q])
                t = (1.0 if theta >= 0 else -1.0) / (abs(theta) + math.sqrt(theta * theta + 1))
                c = 1.0 / math.sqrt(t * t + 1)
                s = t * c
                for k in range(n):
                    akp = a[k][p]
                    akq = a[k][q]
                    a[k][p] = c * akp - s * akq
                    a[k][q] = s * akp + c * akq
                for k in range(n):
                    apk = a[p][k]
                    aqk = a[q][k]
                    a[p][k] = c * apk - s * aqk
                    a[q][k] = s * apk + c * aqk
    return min(a[i][i] for i in range(n))


def gershgorin_lower_bound(matrix):
    """Gershgorin 圆盘下界（对称矩阵特征值 >= min_i(a_ii - R_i)；大矩阵近似检查）。"""
    n = len(matrix)
    bounds = []
    for i in range(n):
        r_i = sum(abs(matrix[i][j]) for j in range(n) if j != i)
        bounds.append(matrix[i][i] - r_i)
    return min(bounds)


def main(argv=None):
    args = parse_args(argv)

    # 顺序文件对账
    order = []
    with open(args.matrix_ids, encoding="utf-8") as handle:
        header = handle.readline()
        if not header.lower().startswith("fid"):
            fail(f"sample_order.tsv 表头异常: {header!r}")
        for line in handle:
            parts = line.split()
            if len(parts) >= 2:
                order.append((parts[0], parts[1]))
    fam_order = []
    with open(args.fam, encoding="utf-8") as handle:
        for line in handle:
            parts = line.split()
            if len(parts) >= 2:
                fam_order.append((parts[0], parts[1]))

    report = {
        "schema": "gwasflow.kinship_validation.v1",
        "matrix": str(Path(args.matrix).resolve()),
        "tolerances": {
            "symmetry_relative": args.symmetry_tolerance,
            "eigenvalue_relative": args.eigenvalue_tolerance,
        },
    }

    # 检查 1: 维度
    matrix = []
    with open(args.matrix, encoding="utf-8") as handle:
        for line in handle:
            row = [float(x) for x in line.split()]
            matrix.append(row)
    n = len(matrix)
    dimension_ok = all(len(row) == n for row in matrix) and n == len(order) == len(fam_order)
    report["n"] = n
    report["dimension_ok"] = dimension_ok
    if not dimension_ok:
        report["failures"] = [
            f"维度不匹配: matrix {n}×{len(matrix[0]) if matrix else 0}, "
            f"order {len(order)}, fam {len(fam_order)}"
        ]
        _write_report_and_fail(args, report)

    # 检查 2: 行列身份（sample_order 与 FAM 顺序逐项一致）
    identity_ok = order == fam_order
    report["identity_ok"] = identity_ok
    if not identity_ok:
        report["failures"] = ["sample_order.tsv 与 FAM 顺序不一致（逐项对账失败）"]
        _write_report_and_fail(args, report)

    # 检查 3: 数值有限性
    finite_ok = all(
        x == x and abs(x) != float("inf")
        for row in matrix for x in row
    )
    report["finite_ok"] = finite_ok
    if not finite_ok:
        report["failures"] = ["矩阵含 NaN/Inf 非有限值（请检查标记集与样本规模）"]
        _write_report_and_fail(args, report)

    # 检查 4: 对称性（容差按矩阵尺度：tol × max(1, |K|max)）
    max_abs = max(abs(x) for row in matrix for x in row)
    scale = max(1.0, max_abs)
    sym_tol_abs = args.symmetry_tolerance * scale
    max_asym = max(abs(matrix[i][j] - matrix[j][i])
                   for i in range(n) for j in range(i + 1, n))
    symmetric_ok = max_asym <= sym_tol_abs
    report["symmetric_ok"] = symmetric_ok
    report["symmetry"] = {
        "max_abs_asymmetry": max_asym,
        "absolute_tolerance": sym_tol_abs,
        "matrix_scale": scale,
    }
    if not symmetric_ok:
        report["failures"] = [
            f"对称性失败: max|K_ij-K_ji|={max_asym} > 容差 {sym_tol_abs}"
            f"（按矩阵尺度 {scale} 声明）"
        ]
        _write_report_and_fail(args, report)

    # 检查 5: PSD（最小特征值；仅容忍数值误差级别，禁止静默修复）
    min_eig = None
    method = None
    try:
        import numpy as np  # noqa: F401

        import numpy
        eigvals = numpy.linalg.eigvalsh(numpy.array(matrix))
        min_eig = float(eigvals[0])
        max_eig = float(eigvals[-1])
        method = "numpy.eigvalsh"
    except ImportError:
        if n <= 800:
            min_eig = jacobi_min_eigenvalue(matrix)
            max_eig = max(matrix[i][i] for i in range(n))
            method = "jacobi_rotation（标准库）"
        else:
            min_eig = gershgorin_lower_bound(matrix)
            max_eig = max(max(row) for row in matrix)
            method = "gershgorin_lower_bound（大矩阵近似下界；结论保守）"

    eig_scale = max(1.0, abs(max_eig) if max_eig is not None else 1.0)
    eig_tol_abs = args.eigenvalue_tolerance * eig_scale
    psd_ok = min_eig is not None and min_eig >= -eig_tol_abs
    report["psd_check"] = {
        "method": method,
        "min_eigenvalue": min_eig,
        "max_eigenvalue": max_eig,
        "absolute_tolerance": -eig_tol_abs,
        "ok": psd_ok,
        "note": "仅容忍数值误差级别；不静默强制修复（谱裁剪属显式后续操作）",
    }
    if not psd_ok:
        report["failures"] = [
            f"PSD 检查失败: min_eigenvalue={min_eig} < -{eig_tol_abs}"
            f"（超出数值误差容忍；请检查标记集与样本规模）"
        ]
        _write_report_and_fail(args, report)

    report["status"] = "passed"
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)
    print(
        f"[validate_kinship] 验证通过: N={n}, 对称={symmetric_ok}, "
        f"min_eig={min_eig}（{method}）"
    )
    return EXIT_OK


def _write_report_and_fail(args, report):
    report["status"] = "failed"
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)
    fail("; ".join(report.get("failures", ["验证失败"])))


if __name__ == "__main__":
    sys.exit(main())
