#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  printf 'Usage: %s <workflow-id> <sample-name>\n' "$0" >&2
  exit 2
fi

workflow_id=$1
sample_name=$2

if [[ ! $workflow_id =~ ^[[:alnum:]-]+$ ]]; then
  printf 'Invalid workflow ID: %s\n' "$workflow_id" >&2
  exit 2
fi

if [[ ! $sample_name =~ ^[[:alnum:]][[:alnum:]._-]*$ ]]; then
  printf 'Invalid sample name: %s\n' "$sample_name" >&2
  exit 2
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "$script_dir/../.." && pwd)
workflow_dir="$repo_root/cromwell-executions/MitochondriaPipeline/$workflow_id"
outputs_workflow_dir="$repo_root/tests/cromwell-outputs/MitochondriaPipeline/$workflow_id"
output_dir="$outputs_workflow_dir/$sample_name"
final_vcf="$workflow_dir/call-LiftOverAfterSelf/execution/out/${sample_name}.self.ref.split.selfToRef.final.vcf"
diagnostic_tsv="$workflow_dir/call-LiftoverSelfCoverage/execution/out/${sample_name}_mtanalysis_diagnostic_statistics.tsv"

if [[ ! -f $final_vcf ]]; then
  printf 'Final VCF not found: %s\n' "$final_vcf" >&2
  exit 1
fi

if [[ ! -f $diagnostic_tsv ]]; then
  printf 'Diagnostic TSV not found: %s\n' "$diagnostic_tsv" >&2
  exit 1
fi

if [[ -e $output_dir ]]; then
  printf 'Output directory already exists; refusing to overwrite: %s\n' "$output_dir" >&2
  exit 1
fi

mkdir -p -- "$outputs_workflow_dir"
staging_dir=$(mktemp -d "$outputs_workflow_dir/.${sample_name}.outputs.XXXXXX")
cleanup() {
  if [[ -n ${staging_dir:-} && -d $staging_dir ]]; then
    rm -f -- "$staging_dir"/*
    rmdir -- "$staging_dir"
  fi
}
trap cleanup EXIT

cp -- "$final_vcf" "$staging_dir/"
cp -- "$diagnostic_tsv" "$staging_dir/"
mv -T -- "$staging_dir" "$output_dir"
staging_dir=

printf 'Copied final VCF and diagnostic TSV to:\n%s\n' "$output_dir"
