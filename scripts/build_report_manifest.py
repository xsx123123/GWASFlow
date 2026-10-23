#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""BioReport 渲染产物清单构建（15.report/report_manifest.json）。

扫描报告输出目录，登记 HTML/PDF 交付件；无产物时非零退出（禁止空清单伪装成功）。
"""

import argparse
import json
import sys
from pathlib import Path

EXIT_OK = 0
EXIT_INPUT_ERROR = 2


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="build_report_manifest.py",
        description="BioReport 渲染产物 → report_manifest.json",
    )
    parser.add_argument("--report-output-dir", required=True, help="bioreport -o 输出目录")
    parser.add_argument("--config", default="", help="report_config.json 路径（登记用）")
    parser.add_argument("--out", required=True, help="report_manifest.json 输出路径")
    args = parser.parse_args(argv)

    out_dir = Path(args.report_output_dir)
    if not out_dir.is_dir():
        print(f"[build_report_manifest ERROR] 报告输出目录不存在: {out_dir}", file=sys.stderr)
        return EXIT_INPUT_ERROR

    artifacts = []
    for path in sorted(out_dir.rglob("*")):
        if path.is_file() and path.suffix in (".html", ".pdf"):
            artifacts.append({
                "path": str(path),
                "format": path.suffix.lstrip("."),
                "size_bytes": path.stat().st_size,
            })
    if not artifacts:
        print(f"[build_report_manifest ERROR] 输出目录无 HTML/PDF 产物: {out_dir}", file=sys.stderr)
        return EXIT_INPUT_ERROR

    manifest = {
        "schema": "gwasflow.report.v1",
        "engine": "bioreport",
        "config": args.config,
        "language_editions": ["zh", "en"],
        "formats": sorted({item["format"] for item in artifacts}),
        "artifacts": artifacts,
        "n_artifacts": len(artifacts),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[build_report_manifest] {len(artifacts)} artifacts -> {out}")
    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
