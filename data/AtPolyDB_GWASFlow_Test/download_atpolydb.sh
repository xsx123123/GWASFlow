#!/usr/bin/env bash
set -euo pipefail

raw_dir="${1:-data/raw}"
mkdir -p "$raw_dir"
base_url="https://raw.githubusercontent.com/Gregor-Mendel-Institute/atpolydb/master"

fetch() {
  local relative_path="$1" file_name="$2" expected_sha="$3"
  if [[ ! -f "$raw_dir/$file_name" ]]; then
    curl --fail --location --retry 3 --connect-timeout 20 \
      "$base_url/$relative_path" --output "$raw_dir/$file_name.part"
    mv "$raw_dir/$file_name.part" "$raw_dir/$file_name"
  fi
  python3 - "$raw_dir/$file_name" "$expected_sha" <<'PY'
import hashlib, pathlib, sys
path = pathlib.Path(sys.argv[1])
digest = hashlib.sha256()
with path.open('rb') as stream:
    for chunk in iter(lambda: stream.read(1024 * 1024), b''):
        digest.update(chunk)
actual = digest.hexdigest()
if actual != sys.argv[2]:
    raise SystemExit(f'SHA256 mismatch: {path}\nExpected {sys.argv[2]}\nObserved {actual}\nKeep the file for investigation; do not bypass this check.')
print(f'OK {path.name} {path.stat().st_size} bytes SHA256={actual}')
PY
}

fetch '250k_snp_data/call_method_75.tar.gz' 'call_method_75.tar.gz' \
  '7ad5ded6c72b703603f3d29ce8cd231b67979dcc56afb462a7035257e752cb76'
fetch 'miscellaneous_data/phenotype_published_raw.tsv' 'phenotype_published_raw.tsv' \
  '6d76c6b59df7b954c1533a383e2322226d684d58919c09020f08e6c0b9004598'

printf '%s\n' 'Downloads verified. Use prepare_atpolydb.py to select TAIR9 and match ecotype IDs.'
