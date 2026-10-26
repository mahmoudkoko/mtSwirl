#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "$script_dir/../.." && pwd)
hail_conda_prefix="${HAIL_CONDA_PREFIX:-$HOME/miniconda3/envs/hail}"
mtswirl_conda_prefix="${MTSWIRL_CONDA_PREFIX:-$HOME/miniconda3/envs/mtSwirl}"
r_conda_prefix="${R_CONDA_PREFIX:-$HOME/miniconda3/envs/R}"
cromwell_jar="$repo_root/resources/cromwell/cromwell-92.jar"
config="$repo_root/inputs/config/cromwell-local-single.conf"
workflow="$repo_root/inputs/wdl/MitoPipeline.wdl"
inputs="$repo_root/inputs/json/prepopulated_inputs_offline_single.json"
copy_script="$script_dir/copy_single_sample_outputs.sh"

if [[ ! -x "$hail_conda_prefix/bin/python" ]]; then
  printf 'Hail Python was not found: %s\n' "$hail_conda_prefix/bin/python" >&2
  printf 'Set HAIL_CONDA_PREFIX to the Hail Conda environment.\n' >&2
  exit 1
fi

if [[ ! -d "$mtswirl_conda_prefix/bin" ]]; then
  printf 'mtSwirl Conda bin directory was not found: %s\n' "$mtswirl_conda_prefix/bin" >&2
  printf 'Set MTSWIRL_CONDA_PREFIX to the mtSwirl Conda environment.\n' >&2
  exit 1
fi

if [[ ! -x "$r_conda_prefix/bin/Rscript" ]]; then
  printf 'Rscript was not found: %s\n' "$r_conda_prefix/bin/Rscript" >&2
  printf 'Set R_CONDA_PREFIX to the R Conda environment.\n' >&2
  exit 1
fi

export PATH="$mtswirl_conda_prefix/bin:$hail_conda_prefix/bin:$r_conda_prefix/bin:$PATH"
export JAVA_HOME="$mtswirl_conda_prefix"
export HAIL_CONDA_PREFIX="$hail_conda_prefix"
export MTSWIRL_CONDA_PREFIX="$mtswirl_conda_prefix"
export R_CONDA_PREFIX="$r_conda_prefix"
export PYSPARK_PYTHON="$hail_conda_prefix/bin/python"

if ! command -v java >/dev/null 2>&1; then
  printf 'Java was not found on PATH.\n' >&2
  exit 1
fi

for file in "$cromwell_jar" "$config" "$workflow" "$inputs" "$copy_script"; do
  if [[ ! -f $file ]]; then
    printf 'Required file not found: %s\n' "$file" >&2
    exit 1
  fi
done

if ! command -v python3 >/dev/null 2>&1; then
  printf 'python3 was not found on PATH; it is required to read workflow metadata and inputs.\n' >&2
  exit 1
fi

metadata_dir=$(mktemp -d "${TMPDIR:-/tmp}/mtSwirl-cromwell.XXXXXX")
metadata_file="$metadata_dir/metadata.json"
cleanup() {
  if [[ -f $metadata_file ]]; then
    rm -f -- "$metadata_file"
  fi
  if [[ -d $metadata_dir ]]; then
    rmdir -- "$metadata_dir"
  fi
}
trap cleanup EXIT

cd -- "$repo_root"
java "-Dconfig.file=$config" -jar "$cromwell_jar" run \
  -m "$metadata_file" \
  "$workflow" \
  -i "$inputs"

readarray -t run_details < <(python3 - "$metadata_file" "$inputs" <<'PY'
import json
import sys

metadata_path, inputs_path = sys.argv[1:]
with open(metadata_path, encoding="utf-8") as metadata_file:
    metadata = json.load(metadata_file)
with open(inputs_path, encoding="utf-8") as inputs_file:
    inputs = json.load(inputs_file)

workflow_id = metadata.get("id")
sample_name = inputs.get("MitochondriaPipeline.sample_name")
if not isinstance(workflow_id, str) or not workflow_id:
    raise SystemExit("Cromwell metadata did not contain a workflow id.")
if not isinstance(sample_name, str) or not sample_name:
    raise SystemExit("The input JSON did not contain a sample name.")

print(workflow_id)
print(sample_name)
PY
)

if [[ ${#run_details[@]} -ne 2 ]]; then
  printf 'Could not read workflow ID and sample name from Cromwell metadata.\n' >&2
  exit 1
fi

bash "$copy_script" "${run_details[0]}" "${run_details[1]}"
