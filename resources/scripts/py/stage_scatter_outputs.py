#!/usr/bin/env python3

import argparse
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path


PREFIX = "SwirlPipeline."


def parse_args():
    parser = argparse.ArgumentParser(description="Stage curated outputs from a successful scatter run.")
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--sample-list", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    return parser.parse_args()


def main():
    args = parse_args()
    metadata = json.loads(args.metadata.read_text(encoding="utf-8"))
    workflow_id = metadata.get("id")
    if not isinstance(workflow_id, str) or not workflow_id:
        raise ValueError("Cromwell metadata did not contain a workflow ID.")
    if metadata.get("status") != "Succeeded":
        raise ValueError(f"Expected a successful Cromwell workflow; got status {metadata.get('status')!r}.")
    outputs = metadata.get("outputs")
    if not isinstance(outputs, dict):
        raise ValueError("Cromwell metadata did not contain workflow outputs.")

    sample_names = [
        line.strip()
        for line in args.sample_list.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    run_dir = args.output_root / workflow_id
    run_dir.parent.mkdir(parents=True, exist_ok=True)
    if run_dir.exists():
        raise FileExistsError(f"Refusing to overwrite existing output directory: {run_dir}")
    staging = Path(tempfile.mkdtemp(prefix=f".{workflow_id}.", dir=run_dir.parent))

    def output_path(name):
        value = outputs.get(PREFIX + name)
        if not isinstance(value, str) or not value:
            raise ValueError(f"Missing Cromwell output {PREFIX}{name}.")
        path = Path(value)
        if not path.is_file():
            raise FileNotFoundError(f"Cromwell output file is unavailable: {path}")
        return path

    def copy_single(name, destination):
        shutil.copy2(output_path(name), staging / destination)

    def copy_per_sample(name, suffix):
        values = outputs.get(PREFIX + name)
        if not isinstance(values, list) or len(values) != len(sample_names):
            raise ValueError(f"Expected {len(sample_names)} files in Cromwell output {name}.")
        for sample, value in zip(sample_names, values):
            if not isinstance(value, str) or not value:
                raise ValueError(f"Missing {name} output for sample {sample}.")
            source = Path(value)
            if not source.is_file():
                raise FileNotFoundError(f"Cromwell output file is unavailable: {source}")
            shutil.copy2(source, staging / f"{sample}.{suffix}")

    try:
        copy_single("merged_calls", "merged_calls.vcf.bgz")
        copy_single("merged_coverage", "merged_coverage.tsv.bgz")
        copy_single("merged_statistics", "merged_statistics.tsv")
        copy_per_sample("final_vcf", "final.vcf")
        copy_per_sample("final_base_level_coverage_metrics", "coverage.tsv")
        copy_per_sample("stats_outputs", "statistics.tsv")
        os.replace(staging, run_dir)
    finally:
        if staging.exists():
            shutil.rmtree(staging)

    print(f"Copied merged outputs and per-sample VCF, coverage, and statistics to:\n{run_dir}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, json.JSONDecodeError) as error:
        print(error, file=sys.stderr)
        raise SystemExit(1) from error
