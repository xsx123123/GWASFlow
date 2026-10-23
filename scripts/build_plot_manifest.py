#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GWASFlow 图表清单构建（plot_manifest.json：与实际图文件一一对应）。"""
import argparse
import json
import os
import sys

EXIT_OK = 0
EXIT_INPUT_ERROR = 2


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="build_plot_manifest.py",
        description="plot_manifest.json：图身份、来源数据、格式（一一对应）。",
    )
    parser.add_argument("--plot-file", action="append", required=True,
                        help="图/绘图数据文件路径（可多次）")
    parser.add_argument("--not-generated", action="append", default=[],
                        help="未生成图及原因（可多次）")
    parser.add_argument("--not-generated-file", default=None,
                        help="未生成图及原因文件（每行一条；避免 shell 引号问题）")
    parser.add_argument("--out", required=True, help="manifest JSON 输出路径")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if args.not_generated_file and os.path.exists(args.not_generated_file):
        with open(args.not_generated_file, encoding="utf-8") as handle:
            args.not_generated.extend(
                line.strip() for line in handle if line.strip()
            )
    entries = []
    missing = []
    for path in args.plot_file:
        if not os.path.exists(path):
            missing.append(path)
            continue
        role = (
            "gwas_plot" if ("manhattan" in path or "/qq." in path)
            else "diagnostic_plot" if ("diagnostics" in path
                                       or "distribution" in path
                                       or "sample_inclusion" in path)
            else "plot_data"
        )
        entries.append({
            "path": path,
            "role": role,
            "format": os.path.splitext(path)[1].lstrip("."),
            "size_bytes": os.path.getsize(path),
        })
    if missing:
        print(f"[build_plot_manifest ERROR] 声明的图文件缺失: {missing}",
              file=sys.stderr)
        sys.exit(EXIT_INPUT_ERROR)
    manifest = {
        "schema": "gwasflow.plot_manifest.v1",
        "plots": entries,
        "n_plots": len(entries),
        "not_generated": list(args.not_generated),
        "note": "图与实际文件一一对应；模块关闭时对应图不生成并标注原因",
    }
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)
    print(f"[build_plot_manifest] {len(entries)} 图登记")
    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
