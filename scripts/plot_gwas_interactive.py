#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GWASFlow 交互式曼哈顿/QQ 图（plotly WebGL，浏览器可交互 HTML）。

实现参照 plotly.com/python/manhattan-plot 教程的交互理念（WebGL 大数据量、
悬停信息、阈值线），因流程契约（p_underflow 显示上限、阈值以 metadata 为
唯一来源、染色体顺序来自 chromosome_table 配置）改用 go.Scattergl 直接构建，
不引入 dash-bio 依赖。

与静态图（plot_gwas.R，svg/pdf/png，发表用）的语义一致性：
- 阈值从 metadata 读取（唯一来源，不另抄）；
- p_underflow 显示上限（p_value_text=="0" 时显示 cap，统计表保留原值）；
- 染色体顺序/长度来自 chromosome_table；表外染色体剔除并审计计数；
- 显著点永不因抽样丢失（交互图默认全量渲染，仅超过 --max-render-points
  时才抽样：显著点全保留 + 非显著均匀随机抽样，seed 可复现）。

交互能力：缩放/平移/框选/套索/悬停查看位点身份与统计量/工具栏导出 PNG。
HTML 默认内嵌 plotly.js（离线可打开；--plotlyjs cdn 可减 ~3.5MB/文件）。
"""
import argparse
import csv
import json
import random
import sys

EXIT_OK = 0
EXIT_INPUT_ERROR = 2

# 与静态图一致的配色（chromosome 交替 + 显著/阈值）
CHROM_COLORS = ["#2C6E9B", "#D1825E"]
SIG_COLOR = "#C0392B"
NSIG_COLOR = "#2C6E9B"
THRESH_COLOR = "#C0392B"
CAP_COLOR = "#7F8C8D"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="plot_gwas_interactive.py",
        description="GWASFlow 交互式曼哈顿/QQ 图（plotly Scattergl → HTML）",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--sumstats", required=True, help="08.results/{trait}/gwas.tsv")
    parser.add_argument("--metadata", required=True,
                        help="08.results/{trait}/metadata.json（阈值唯一来源）")
    parser.add_argument("--chromosome-table", default=None,
                        help="染色体顺序+长度表（两列带表头：名称、长度 bp）")
    parser.add_argument("--trait", required=True)
    parser.add_argument("--model", default="gemma_lmm")
    parser.add_argument("--max-render-points", type=int, default=2000000,
                        help="交互图全量渲染上限；超过才抽样（显著点全保留）")
    parser.add_argument("--p-zero-display-cap", default="1e-300",
                        help="P=0 下溢显示上限（统计表保留原值，不修改）")
    parser.add_argument("--seed", type=float, default=20261023, help="抽样种子（可复现）")
    parser.add_argument("--plotlyjs", choices=["embed", "cdn"], default="embed",
                        help="HTML 内嵌 plotly.js（离线可开）或 CDN 引用（文件更小）")
    parser.add_argument("--out-dir", required=True, help="输出目录")
    return parser.parse_args(argv)


def fail(message, code=EXIT_INPUT_ERROR):
    print(f"[plot_gwas_interactive ERROR] {message}", file=sys.stderr)
    sys.exit(code)


def main(argv=None):
    args = parse_args(argv)

    try:
        import plotly.graph_objects as go
    except ImportError:
        fail("plotly 未安装（envs/gwasflow-plotly.yaml；或 pip install plotly）")

    # ---- metadata（阈值唯一来源）----
    with open(args.metadata, encoding="utf-8") as handle:
        meta = json.load(handle)
    threshold_text = str(meta.get("multiple_testing", {}).get("threshold", "")).strip()
    m_tests = meta.get("counts", {}).get("n_valid_p")
    try:
        threshold_num = float(threshold_text)
    except ValueError:
        threshold_num = None
        print(f"[plot_gwas_interactive WARNING] metadata 阈值不可解析: {threshold_text!r}，不画阈值线")

    # ---- 染色体表（顺序 + 长度）----
    chrom_order, length_map = [], {}
    if args.chromosome_table:
        with open(args.chromosome_table, newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            cols = reader.fieldnames or []
            if len(cols) < 2:
                fail(f"chromosome_table 列数不足: {cols}")
            for row in reader:
                name = str(row[cols[0]]).strip()
                chrom_order.append(name)
                try:
                    length_map[name] = float(row[cols[1]])
                except (TypeError, ValueError):
                    fail(f"chromosome_table 长度非数值: {row[cols[1]]!r} @ {name}")

    # ---- 统计表加载 ----
    rows = []
    with open(args.sumstats, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for col in ("chrom", "pos", "variant_id", "p_value", "p_value_text",
                    "beta", "effect_allele", "significant"):
            if col not in (reader.fieldnames or []):
                fail(f"gwas.tsv 缺失列: {col}")
        for row in reader:
            rows.append(row)

    # ---- 显示用 -log10p（与静态图同语义；不改统计表原始值）----
    try:
        p_cap = float(args.p_zero_display_cap)
    except ValueError:
        fail(f"--p-zero-display-cap 非数值: {args.p_zero_display_cap!r}")
    import math
    if p_cap <= 0:
        fail(f"--p-zero-display-cap 必须为正数: {args.p_zero_display_cap!r}")
    cap_neglog = -math.log10(p_cap)

    records = []
    n_outside_table = 0
    for row in rows:
        chrom = str(row["chrom"]).strip()
        if chrom_order and chrom not in chrom_order:
            n_outside_table += 1
            continue
        try:
            p_num = float(row["p_value"])
        except (TypeError, ValueError):
            p_num = None
        is_underflow = (row.get("p_value_text", "") or "").strip() == "0"
        if p_num is not None and p_num > 0:
            neglog = -math.log10(p_num)
        else:
            neglog = None
        neglog_display = cap_neglog if is_underflow else neglog
        significant = (row.get("significant", "") or "").strip().lower() == "true"
        if not significant and threshold_num is not None and p_num is not None:
            significant = p_num > 0 and p_num <= threshold_num
        try:
            beta = float(row["beta"])
            beta_text = f"{beta:.6g}"
        except (TypeError, ValueError):
            beta_text = "NA"
        p_hover = (f"0 (< display cap {args.p_zero_display_cap})"
                   if is_underflow else (row.get("p_value_text") or row["p_value"]))
        records.append({
            "chrom": chrom, "pos": int(row["pos"]), "vid": row["variant_id"],
            "neglog_display": neglog_display, "neglog": neglog,
            "significant": significant, "underflow": is_underflow,
            "beta_text": beta_text, "eff": row["effect_allele"], "p_hover": p_hover,
        })

    n_total = len(records)
    n_underflow = sum(1 for r in records if r["underflow"])
    n_significant = sum(1 for r in records if r["significant"])

    # ---- 抽样（仅超 --max-render-points；显著点全保留，seed 可复现）----
    sampled_note = None
    if n_total > args.max_render_points:
        rng = random.Random(args.seed)
        sig_idx = [i for i, r in enumerate(records) if r["significant"]]
        non_sig_idx = [i for i, r in enumerate(records)
                       if not r["significant"] and r["neglog_display"] is not None]
        budget = max(0, args.max_render_points - len(sig_idx))
        keep = set(sig_idx) | set(rng.sample(non_sig_idx, min(budget, len(non_sig_idx))))
        keep |= {i for i, r in enumerate(records) if r["neglog_display"] is None}
        records = [records[i] for i in sorted(keep)]
        sampled_note = (f"非显著抽样 {len(keep) - len(sig_idx)}/{len(non_sig_idx)}"
                        f"（渲染上限 {args.max_render_points}，seed {args.seed}；显著点全保留）")

    # ---- 基因组累计坐标（chromosome_table 长度；与静态图一致）----
    if not chrom_order:
        # 无表时退化为按出现顺序的去重染色体（长度用最大 pos 估计，仅影响显示）
        seen = []
        for r in records:
            if r["chrom"] not in seen:
                seen.append(r["chrom"])
        chrom_order = seen
        length_map = {c: max((r["pos"] for r in records if r["chrom"] == c),
                             default=1) for c in seen}
        print("[plot_gwas_interactive WARNING] 未提供 chromosome_table，"
              "染色体顺序与长度由数据推断（仅显示用途）")
    offsets = [0.0]
    for c in chrom_order:
        offsets.append(offsets[-1] + length_map.get(c, 0))
    chrom_mid = [(offsets[i] + offsets[i + 1]) / 2 for i in range(len(chrom_order))]
    off_by_chrom = dict(zip(chrom_order, offsets[:-1]))

    plot_recs = [r for r in records if r["neglog_display"] is not None]
    n_plot = len(plot_recs)

    # ===================== Manhattan（交互）=====================
    hover = ("<b>%{customdata[0]}</b><br>chr%{customdata[1]}:%{customdata[2]:,d}"
             "<br>P = %{customdata[3]}"
             "<br>beta = %{customdata[4]} (%{customdata[5]} allele)"
             "<extra></extra>")
    fig_m = go.Figure()
    for ci, chrom in enumerate(chrom_order):
        pts = [r for r in plot_recs if r["chrom"] == chrom]
        if not pts:
            continue
        fig_m.add_trace(go.Scattergl(
            x=[off_by_chrom[chrom] + r["pos"] for r in pts],
            y=[r["neglog_display"] for r in pts],
            customdata=[[r["vid"], r["chrom"], r["pos"], r["p_hover"],
                         r["beta_text"], r["eff"]] for r in pts],
            mode="markers",
            marker={"size": 4, "color": CHROM_COLORS[ci % 2],
                    "line": {"width": 0}},
            name=f"chr{chrom}",
            hovertemplate=hover,
        ))

    title_m = (f"Manhattan plot (interactive) — trait: {args.trait}"
               f" (model: {args.model})")
    sub_m = (f"Bonferroni threshold = {threshold_text}"
             f" (alpha/(family_size × M), M = {m_tests})"
             f" | {n_plot} points rendered"
             + (f" | {sampled_note}" if sampled_note else "")
             + (f" | P=0 display cap = {args.p_zero_display_cap}（统计表保留原值）"
                if n_underflow else ""))
    fig_m.update_layout(
        title={"text": f"{title_m}<br><sup>{sub_m}</sup>", "x": 0.5},
        xaxis={"title": "Chromosome",
               "tickvals": chrom_mid, "ticktext": [f"chr{c}" for c in chrom_order],
               "showgrid": False},
        yaxis={"title": "-log10(P)", "showgrid": True, "gridcolor": "#EEEEEE"},
        template="plotly_white",
        hovermode="closest",
        showlegend=False,
        dragmode="zoom",
    )
    if threshold_num is not None and threshold_num > 0:
        fig_m.add_hline(y=-math.log10(threshold_num), line_dash="dash",
                        line_color=THRESH_COLOR, line_width=1,
                        annotation_text=f"Bonferroni {threshold_num:.3g}",
                        annotation_position="top left")
    if n_underflow:
        fig_m.add_hline(y=cap_neglog, line_dash="dot",
                        line_color=CAP_COLOR, line_width=1,
                        annotation_text=f"P=0 display cap ({args.p_zero_display_cap})",
                        annotation_position="bottom left")

    # ===================== QQ（交互）=====================
    qq_recs = [r for r in records if r["neglog"] is not None]
    qq_recs.sort(key=lambda r: (r["neglog"], r["vid"]))
    n_qq = len(qq_recs)
    expected = [-math.log10((i + 1) / (n_qq + 1)) for i in range(n_qq)]
    fig_q = go.Figure()
    for is_sig, color, name in ((False, NSIG_COLOR, "not significant"),
                                (True, SIG_COLOR, "significant")):
        pts = [(e, r) for e, r in zip(expected, qq_recs) if r["significant"] == is_sig]
        if not pts:
            continue
        fig_q.add_trace(go.Scattergl(
            x=[e for e, _ in pts],
            y=[r["neglog_display"] for _, r in pts],
            customdata=[[r["vid"], r["p_hover"], r["chrom"], r["pos"]]
                        for _, r in pts],
            mode="markers",
            marker={"size": 4, "color": color, "line": {"width": 0}},
            name=name,
            hovertemplate=("<b>%{customdata[0]}</b><br>chr%{customdata[2]}"
                           ":%{customdata[3]:,d}<br>P = %{customdata[1]}"
                           "<br>expected -log10(P) = %{x:.3g}"
                           "<br>observed -log10(P) = %{y:.3g}<extra></extra>"),
        ))
    axis_max = max([r["neglog_display"] for r in qq_recs] + expected + [1.0])
    fig_q.update_layout(
        title={"text": (f"QQ plot (interactive) — trait: {args.trait}"
                        f"<br><sup>{n_qq} points（expected = -log10(rank/(n+1))"
                        f"{'；' + sampled_note if sampled_note else ''}）</sup>"),
               "x": 0.5},
        xaxis={"title": "Expected -log10(P)", "range": [0, axis_max * 1.05],
               "showgrid": True, "gridcolor": "#EEEEEE"},
        yaxis={"title": "Observed -log10(P)", "range": [0, axis_max * 1.05],
               "showgrid": True, "gridcolor": "#EEEEEE"},
        template="plotly_white",
        hovermode="closest",
    )
    fig_q.add_shape(type="line", x0=0, y0=0, x1=axis_max * 1.05, y1=axis_max * 1.05,
                    line={"dash": "dash", "color": THRESH_COLOR, "width": 1})

    # ---- 输出 ----
    import os
    os.makedirs(args.out_dir, exist_ok=True)
    plotlyjs = True if args.plotlyjs == "embed" else "cdn"
    html_config = {
        "displaylogo": False,
        "toImageButtonOptions": {"format": "png", "scale": 2},
    }
    manhattan_path = os.path.join(args.out_dir, "manhattan_interactive.html")
    qq_path = os.path.join(args.out_dir, "qq_interactive.html")
    fig_m.write_html(manhattan_path, include_plotlyjs=plotlyjs,
                     config=html_config, full_html=True)
    fig_q.write_html(qq_path, include_plotlyjs=plotlyjs,
                     config=html_config, full_html=True)

    print(
        f"[plot_gwas_interactive] {args.trait}: Manhattan/QQ 交互图生成"
        f"（{n_plot} 点渲染，显著 {n_significant} 全保留"
        + (f"，{sampled_note}" if sampled_note else "，未触发抽样")
        + f"，P=0 下溢 {n_underflow} 个"
        + (f"，表外染色体剔除 {n_outside_table} 行" if n_outside_table else "")
        + f"）→ {args.out_dir}"
    )
    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
