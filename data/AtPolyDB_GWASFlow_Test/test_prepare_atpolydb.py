#!/usr/bin/env python3
"""Contract tests for ecotype matching, diploid calls, assembly choice and failures."""
import io
import json
import subprocess
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).with_name("prepare_atpolydb.py")


class PrepareTests(unittest.TestCase):
    def fixture(self, root, member="call_method_75/call_method_75_TAIR9.csv", bad_call=False, duplicate=False):
        ids = ["104", "101", "102", "103"]
        if duplicate:
            ids[-1] = "104"
        lines = ["-,-,1,2,3,4", "Chromosome,Positions," + ",".join(ids)]
        for chrom in range(1, 6):
            # The NA-phenotype sample 102 carries the only G: this site becomes monomorphic.
            lines.append(f"{chrom},100,A,A,G,A")
            calls = "R,C,T,A" if bad_call else "A,C,T,A"
            lines.append(f"{chrom},200,{calls}")
        data = ("\n".join(lines) + "\n").encode()
        archive = root / "source.tar.gz"
        with tarfile.open(archive, "w:gz") as tf:
            info = tarfile.TarInfo(member)
            info.size = len(data)
            tf.addfile(info, io.BytesIO(data))
        phen = root / "phenotype.tsv"
        phen.write_bytes("ecotype id\tnativename\t5_FT10\n103\tBör\t20\n104\tX\t30\n102\tY\tNA\n101\tZ\t40\n999\tQ\t50\n".encode("latin-1"))
        return archive, phen

    def run_prepare(self, root, archive, phen):
        return subprocess.run([sys.executable, str(SCRIPT), "--archive", str(archive),
                               "--phenotypes", str(phen), "--min-samples", "3",
                               "--out-dir", str(root / "out")], capture_output=True, text=True)

    def test_matching_and_homozygous_diploid_calls(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = self.run_prepare(root, *self.fixture(root))
            self.assertEqual(result.returncode, 0, result.stderr)
            out = root / "out"
            self.assertEqual((out / "samples.tsv").read_text().splitlines()[1:],
                             ["104\t0\t104\ttrue", "101\t0\t101\ttrue", "103\t0\t103\ttrue"])
            self.assertEqual((out / "phenotypes.tsv").read_text(), "sample\tFT10\n104\t30\n101\t40\n103\t20\n")
            tped = [line.split() for line in (out / "atpolydb.tped").read_text().splitlines()]
            self.assertEqual(len(tped), 5)
            self.assertEqual(tped[0], ["1", "1:200", "0", "200", "A", "A", "C", "C", "A", "A"])
            report = json.loads((out / "preparation.json").read_text())
            self.assertEqual(report["selected_nonmissing_samples"], 3)
            self.assertEqual(report["snps"]["skipped_monomorphic"], 5)
            self.assertEqual(report["unmatched_phenotype_ids"], ["999"])
            self.assertEqual(report["matched_missing_phenotype_ids"], ["102"])
            self.assertEqual(report["phenotype_encoding"], "latin-1")

    def test_rejects_unknown_genotype_encoding(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = self.run_prepare(root, *self.fixture(root, bad_call=True))
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Unsupported genotype tokens", result.stderr)
            self.assertFalse((root / "out" / "preparation.json").exists())

    def test_rejects_duplicate_ecotype_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = self.run_prepare(root, *self.fixture(root, duplicate=True))
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Repeated genotype ecotype IDs", result.stderr)

    def test_never_substitutes_tair8(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = self.run_prepare(root, *self.fixture(root, member="call_method_75/call_method_75_TAIR8.csv"))
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("never substitute TAIR8", result.stderr)


if __name__ == "__main__":
    unittest.main()
