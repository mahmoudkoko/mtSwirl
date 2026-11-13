#!/usr/bin/env python3
"""Merge per-sample mtDNA outputs into one set of files per batch.

Repeats the steps of the MergeMitoMultiSampleOutputs WDL task
(v2.6_MongoSwirl_Multi/fullMitoPipeline_v2_6_Multi.wdl) for every batch in a manifest:
  - statistics: row-bound into one table (R rbind in the WDL)
  - coverage:   merge_per_batch.py --run-coverage (Hail)
  - VCF calls:  merge_per_batch.py --run-variants (Hail)

Manifest: tab-separated with one row per sample and the header
  batch<TAB>sample<TAB>vcf<TAB>coverage<TAB>stats
Paths may be local or gs:// paths readable by Hail.

Outputs are written to OUTPUT_DIR/BATCH/ using the same file names as the WDL:
  batch_analysis_statistics.tsv, batch_merged_mt_coverage.tsv.bgz, batch_merged_mt_calls.vcf.bgz
"""
import argparse
import os
import subprocess
import sys

import pandas as pd

REQUIRED_COLUMNS = ["batch", "sample", "vcf", "coverage", "stats"]
FILE_COLUMNS = ["vcf", "coverage", "stats"]
DEFAULT_MERGE_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "merge_per_batch.py")


def read_manifest(path):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        sys.exit(f"Manifest {path} is missing columns: {', '.join(missing)}")

    empty = df[(df[REQUIRED_COLUMNS] == "").any(axis=1)]
    if not empty.empty:
        sys.exit(f"Manifest has rows with empty required fields:\n{empty.to_string()}")

    dup = df[df.duplicated(["batch", "sample"], keep=False)]
    if not dup.empty:
        sys.exit(f"Manifest has duplicate batch/sample pairs:\n{dup.to_string()}")
    return df


def check_local_paths(df):
    missing = sorted({p for col in FILE_COLUMNS for p in df[col]
                      if not p.startswith("gs://") and not os.path.exists(p)})
    if missing:
        sys.exit("Missing input files:\n" + "\n".join(missing))


def write_paths_tsv(samples, paths, out_path):
    # Same layout as the WDL: columns s and path, no index
    pd.DataFrame({"s": samples, "path": paths}).to_csv(out_path, sep="\t", index=False)


def merge_statistics(paths, out_path):
    tables = [pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False) for p in paths]
    columns = list(tables[0].columns)
    for path, table in zip(paths, tables):
        if list(table.columns) != columns:
            sys.exit(f"Statistics file {path} has different columns than {paths[0]}")
    pd.concat(tables, ignore_index=True).to_csv(out_path, sep="\t", index=False)


def run_hail_merge(args, mode, input_tsv, output_path, temp_dir):
    cmd = [
        sys.executable, args.merge_script, mode,
        "--input-tsv", input_tsv,
        "--temp-dir", temp_dir,
        "--output-flat-file", output_path,
        "--n-read-partitions", str(args.n_read_partitions),
        "--chunk-size", str(args.chunk_size),
    ]
    subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest", required=True, help="Tab-separated manifest described above.")
    parser.add_argument("--output-dir", required=True, help="Directory where per-batch outputs are written.")
    parser.add_argument("--batches", nargs="+", help="Only merge these batch IDs (default: all batches).")
    parser.add_argument("--merge-script", default=DEFAULT_MERGE_SCRIPT,
                        help="Path to merge_per_batch.py (default: alongside this script).")
    parser.add_argument("--n-read-partitions", type=int, default=1,
                        help="Partitions used when reading coverage files (as in the WDL).")
    parser.add_argument("--chunk-size", type=int, default=100,
                        help="Number of samples merged per Hail job (as in the WDL).")
    args = parser.parse_args()

    manifest = read_manifest(args.manifest)
    if args.batches:
        unknown = sorted(set(args.batches) - set(manifest["batch"]))
        if unknown:
            sys.exit(f"Batches not found in manifest: {', '.join(unknown)}")
        manifest = manifest[manifest["batch"].isin(args.batches)]
    check_local_paths(manifest)

    for batch, rows in manifest.groupby("batch", sort=True):
        out_dir = os.path.join(args.output_dir, batch)
        temp_dir = os.path.join(out_dir, "tmp")
        os.makedirs(temp_dir, exist_ok=True)
        print(f"[{batch}] merging {len(rows)} samples into {out_dir}", flush=True)

        samples = rows["sample"].tolist()
        coverage_tsv = os.path.join(out_dir, "coverage_paths.tsv")
        vcf_tsv = os.path.join(out_dir, "vcf_paths.tsv")
        write_paths_tsv(samples, rows["coverage"].tolist(), coverage_tsv)
        write_paths_tsv(samples, rows["vcf"].tolist(), vcf_tsv)

        merge_statistics(rows["stats"].tolist(), os.path.join(out_dir, "batch_analysis_statistics.tsv"))
        run_hail_merge(args, "--run-coverage", coverage_tsv,
                       os.path.join(out_dir, "batch_merged_mt_coverage.tsv.bgz"), temp_dir)
        run_hail_merge(args, "--run-variants", vcf_tsv,
                       os.path.join(out_dir, "batch_merged_mt_calls.vcf.bgz"), temp_dir)
        print(f"[{batch}] done", flush=True)


if __name__ == "__main__":
    main()
