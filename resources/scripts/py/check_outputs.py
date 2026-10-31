#!/usr/bin/env python3
"""Smoke-check a staged SwirlPipeline output directory.

Verifies the expected files exist and are non-empty, and reports the number of
variant records per sample VCF. Optionally fails if a count is zero.
"""
import argparse
import gzip
import sys
from pathlib import Path

MERGED = ["merged_calls.vcf.bgz", "merged_coverage.tsv.bgz", "merged_statistics.tsv"]


def count_records(path):
    opener = gzip.open if path.suffix in {".gz", ".bgz"} else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--samples", nargs="+", help="Expected sample IDs (default: infer from *.final.vcf)")
    args = parser.parse_args()

    out = args.output_dir
    samples = args.samples or sorted(p.name[: -len(".final.vcf")] for p in out.glob("*.final.vcf"))
    problems = []
    if not samples:
        problems.append("no *.final.vcf files found")

    expected = list(MERGED)
    for sample in samples:
        expected += [f"{sample}.final.vcf", f"{sample}.coverage.tsv", f"{sample}.statistics.tsv"]
    for name in expected:
        path = out / name
        if not path.is_file() or path.stat().st_size == 0:
            problems.append(f"missing or empty: {name}")

    for sample in samples:
        vcf = out / f"{sample}.final.vcf"
        if vcf.is_file():
            records = count_records(vcf)
            print(f"{sample}: {records} VCF records")
            if records == 0:
                problems.append(f"{sample}: final VCF has no records")

    if problems:
        print("FAILED:\n  " + "\n  ".join(problems), file=sys.stderr)
        return 1
    print("OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
