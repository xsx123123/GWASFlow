#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GWASFlow 复制交付（参数手册 §9.9 契约，D6：Python 标准库自研交付器）。

- 按显式 artifact_index 逐文件复制 + 逐文件 MD5（目录型产物计入 manifest，
  禁止静默跳过；禁止递归 glob）
- 默认 mode=copy（交付目录自洽可独立移动，禁止默认 symlink）
- 输出三件套：delivery_manifest.json / delivery_manifest.md5 / delivery_details.log
- manifest 记 source/delivered_path/size/md5/action=copy
- {input:q} 式安全传参（路径空格/argv 长度问题由 manifest 文件解决）

错误码：0=成功；2=格式/配置/输入错误（索引缺失/源文件缺失）。
"""

import argparse
import csv
import hashlib
import json
import os
import shutil
import sys
from pathlib import Path

EXIT_OK = 0
EXIT_INPUT_ERROR = 2


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="deliver_results.py",
        description="按显式索引复制交付（copy 模式；manifest 三件套）。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--artifact-index", required=True,
                        help="13.summary/artifact_index.tsv（显式集合，不用递归 glob）")
    parser.add_argument("--source-dir", required=True, help="源目录（workflow）")
    parser.add_argument("--target-dir", required=True, help="交付目录（data_deliver）")
    parser.add_argument("--mode", default="copy", choices=["copy"],
                        help="交付模式（默认 copy：可独立移动；symlink 被拒绝）")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--out-dir", required=True, help="manifest 三件套输出目录")
    return parser.parse_args(argv)


def fail(message, code=EXIT_INPUT_ERROR):
    print(f"[deliver_results ERROR] {message}", file=sys.stderr)
    sys.exit(code)


def md5_of_file(path):
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main(argv=None):
    args = parse_args(argv)

    index_path = Path(args.artifact_index)
    if not index_path.exists():
        fail(f"索引文件缺失（禁止静默跳过交付）: {index_path}")

    source_dir = Path(args.source_dir).resolve()
    target_dir = Path(args.target_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    target_dir.mkdir(parents=True, exist_ok=True)

    # ---- 读取显式索引 ----
    files = []
    with open(index_path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "path" not in (reader.fieldnames or []):
            fail(f"artifact_index 缺失 path 列（实际列: {reader.fieldnames}）")
        for line_no, row in enumerate(reader, start=2):
            rel = (row.get("path") or "").strip()
            if not rel:
                continue
            files.append({
                "path": rel,
                "role": (row.get("role") or "").strip(),
                "expected_md5": (row.get("md5") or "").strip(),
            })
    if not files:
        fail("artifact_index 为空（无文件可交付）")

    # ---- 逐文件复制 + MD5 ----
    manifest_files = []
    details_lines = []
    failures = []
    for entry in files:
        rel = entry["path"]
        src = source_dir / rel
        if not src.exists():
            failures.append(f"源文件缺失: {src}")
            continue
        dst = target_dir / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        observed_md5 = md5_of_file(dst)
        if entry["expected_md5"] and entry["expected_md5"] != observed_md5:
            failures.append(
                f"复制后 MD5 不一致: {rel}"
                f"（期望 {entry['expected_md5']}，实际 {observed_md5}）"
            )
        manifest_files.append({
            "source": str(src),
            "delivered_path": str(dst),
            "relative_path": rel,
            "role": entry["role"],
            "size_bytes": os.path.getsize(dst),
            "md5": observed_md5,
            "action": "copy",
        })
        details_lines.append(f"{observed_md5}  {rel}")

    if failures:
        for f in failures:
            print(f"[deliver_results ERROR] {f}", file=sys.stderr)
        fail(f"交付失败（{len(failures)} 项）")

    # ---- manifest 三件套 ----
    manifest = {
        "schema": "gwasflow.delivery.v1",
        "mode": args.mode,
        "source_dir": str(source_dir),
        "target_dir": str(target_dir),
        "n_files": len(manifest_files),
        "files": manifest_files,
        "note": "copy 模式：交付目录自洽可独立移动（禁止默认 symlink）",
    }
    manifest_path = out_dir / "delivery_manifest.json"
    with open(manifest_path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)

    manifest_md5 = md5_of_file(manifest_path)
    with open(out_dir / "delivery_manifest.md5", "w", encoding="utf-8") as handle:
        handle.write(f"{manifest_md5}  delivery_manifest.json\n")

    with open(out_dir / "delivery_details.log", "w", encoding="utf-8") as handle:
        handle.write("\n".join(details_lines) + "\n")

    print(
        f"[deliver_results] 交付完成：{len(manifest_files)} 文件（copy 模式），"
        f"目录 {target_dir}"
    )
    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
