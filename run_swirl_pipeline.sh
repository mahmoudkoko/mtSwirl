#!/usr/bin/env bash
set -euo pipefail

usage() {
  printf 'Usage: %s [--validate-only] [--test-run | --sample-manifest FILE] [--max-concurrent-jobs N] [--cleanup-executions] [--configs-json FILE] [--files-json FILE] [--scripts-json FILE]\n' "$0" >&2
  printf 'Use --sample-manifest to specify a CSV, TSV, or TXT sample manifest.\n' >&2
  printf 'Use --test-run to run with the bundled one-sample manifest at inputs/test/one_test_sample.tsv\n  (a two-sample manifest is at inputs/test/two_test_samples.tsv; pass it with --sample-manifest).\n' >&2
  printf 'Use --max-concurrent-jobs N to override the local Cromwell task concurrency limit.\n' >&2
  printf 'Use --cleanup-executions to delete the run'"'"'s Cromwell execution directory after a successful run and staging.\n' >&2
  printf 'Uses the packaged JSON files configured in inputs/json/pipeline_configs.json by default.\n' >&2
}

validate_only=false
test_run=false
configs_json=""
files_json=""
scripts_json=""
manifest=""
manifest_seen=false
max_concurrent_jobs=""
cleanup_executions=false

while (($#)); do
  case "$1" in
    --validate-only)
      validate_only=true
      shift
      ;;
    --cleanup-executions)
      cleanup_executions=true
      shift
      ;;
    --test-run)
      test_run=true
      shift
      ;;
    --sample-manifest)
      if (($# < 2)); then
        printf 'Option --sample-manifest requires a file path.\n' >&2
        usage
        exit 2
      fi
      if [[ $manifest_seen == true ]]; then
        printf 'Only one sample manifest may be specified.\n' >&2
        usage
        exit 2
      fi
      manifest=$2
      manifest_seen=true
      shift 2
      ;;
    --max-concurrent-jobs)
      if (($# < 2)); then
        printf 'Option --max-concurrent-jobs requires a positive integer.\n' >&2
        usage
        exit 2
      fi
      max_concurrent_jobs=$2
      shift 2
      ;;
    --configs-json|--files-json|--scripts-json)
      if (($# < 2)); then
        usage
        exit 2
      fi
      case "$1" in
        --configs-json) configs_json=$2 ;;
        --files-json) files_json=$2 ;;
        --scripts-json) scripts_json=$2 ;;
      esac
      shift 2
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    -*)
      printf 'Unknown option: %s\n' "$1" >&2
      usage
      exit 2
      ;;
    *)
      printf 'Unexpected positional argument: %s\nUse --sample-manifest FILE to specify the manifest.\n' "$1" >&2
      usage
      exit 2
      ;;
  esac
done

if [[ -n $max_concurrent_jobs ]] &&
  { [[ ! $max_concurrent_jobs =~ ^[1-9][0-9]*$ ]] ||
    ((${#max_concurrent_jobs} > 10)) ||
    ((10#$max_concurrent_jobs > 20)); }; then
  printf 'Option --max-concurrent-jobs must be an integer between 1 and 20.\n' >&2
  exit 2
fi

if [[ $test_run == true && $manifest_seen == true ]]; then
  printf 'Do not combine --test-run with --sample-manifest; --test-run uses the bundled test manifest.\n' >&2
  usage
  exit 2
fi

if [[ $test_run == false && $manifest_seen == false ]]; then
  usage
  exit 2
fi

script_path=$(readlink -f -- "${BASH_SOURCE[0]}")
repo_root=$(cd -- "$(dirname -- "$script_path")" && pwd)
if [[ $test_run == true ]]; then
  manifest="$repo_root/inputs/test/one_test_sample.tsv"
fi
manifest=$(realpath -- "$manifest")
configs_json=${configs_json:-"$repo_root/inputs/json/pipeline_configs.json"}
configs_json=$(realpath -- "$configs_json")
config_reader="$repo_root/resources/scripts/py/read_scatter_config.py"

for file in "$manifest" "$configs_json" "$config_reader"; do
  if [[ ! -f $file ]]; then
    printf 'Required input file not found: %s\n' "$file" >&2
    exit 1
  fi
done

if ! command -v python3 >/dev/null 2>&1; then
  printf 'python3 was not found on PATH; it is required to read pipeline configuration.\n' >&2
  exit 1
fi

config_args=(--config-json "$configs_json" --repo-root "$repo_root")
if [[ -n $files_json ]]; then
  config_args+=(--files-json "$files_json")
fi
if [[ -n $scripts_json ]]; then
  config_args+=(--scripts-json "$scripts_json")
fi
config_output=$(python3 "$config_reader" "${config_args[@]}")
mapfile -t config_values <<< "$config_output"
if ((${#config_values[@]} != 8)); then
  printf 'Pipeline config did not return all 8 expected values.\n' >&2
  exit 1
fi

tools_path=${config_values[0]}
workflow=${config_values[1]}
prepare_script=${config_values[2]}
stage_script=${config_values[3]}
output_root=${config_values[4]}
files_json=${config_values[5]}
scripts_json=${config_values[6]}
haplocheck_cli=${config_values[7]}

for file in "$workflow" "$prepare_script" "$stage_script" "$files_json" "$scripts_json" "$haplocheck_cli"; do
  if [[ ! -f $file ]]; then
    printf 'Configured file not found: %s\n' "$file" >&2
    exit 1
  fi
done
if [[ -n $tools_path ]]; then
  for bin_dir in "$tools_path" "$tools_path/bin"; do
    [[ -d $bin_dir ]] && PATH="$bin_dir:$PATH"
  done
  export PATH
fi

for tool in java Rscript cromwell womtool; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    printf '%s was not found in tools_path or PATH.\n' "$tool" >&2
    exit 1
  fi
done

if ! python -c 'import hail' >/dev/null 2>&1; then
  printf 'python on PATH (%s) cannot import hail; activate or configure the environment with Hail.\n' "$(command -v python || echo none)" >&2
  exit 1
fi

java_home_for() {
  local resolved_java java_dir
  resolved_java=$(readlink -f -- "$1")
  java_dir=$(dirname -- "$resolved_java")
  if [[ $(basename -- "$java_dir") == bin ]]; then
    dirname -- "$java_dir"
  else
    printf '%s\n' "$java_dir"
  fi
}

export JAVA_HOME
JAVA_HOME=$(java_home_for "$(command -v java)")
export PYSPARK_PYTHON
PYSPARK_PYTHON=$(command -v python)

mkdir -p -- "$repo_root/outputs/tmp"
metadata_dir=$(mktemp -d "$repo_root/outputs/tmp/mtSwirl-scatter.XXXXXX")
metadata_file="$metadata_dir/metadata.json"
inputs_file="$metadata_dir/inputs.json"
cram_list="$metadata_dir/crams.list"
index_list="$metadata_dir/indexes.list"
sample_list="$metadata_dir/samples.list"
cleanup() {
  rm -f -- "$metadata_file" "$inputs_file" "$cram_list" "$index_list" "$sample_list"
  rmdir -- "$metadata_dir"
  rmdir -- "$repo_root/outputs/tmp" 2>/dev/null || true
}
trap cleanup EXIT

python3 "$prepare_script" \
  --manifest "$manifest" \
  --repo-root "$repo_root" \
  --files-json "$files_json" \
  --scripts-json "$scripts_json" \
  --configs-json "$configs_json" \
  --haplocheck-cli "$haplocheck_cli" \
  --cram-list "$cram_list" \
  --index-list "$index_list" \
  --sample-list "$sample_list" \
  --inputs-output "$inputs_file"

cd -- "$repo_root"
womtool validate -i "$inputs_file" -l "$workflow"
if [[ $validate_only == true ]]; then
  printf 'Validation only; Cromwell was not started.\n'
  exit 0
fi

cromwell_args=(
  "-Dworkflow-options.workflow-log-dir=$repo_root/outputs/cromwell-logs"
  "-Dworkflow-options.workflow-log-temporary=false"
  "-Dbackend.providers.Local.config.root=$repo_root/outputs/cromwell-executions"
  "-Dbackend.providers.Local.config.concurrent-job-limit=${max_concurrent_jobs:-1}"
)
cromwell "${cromwell_args[@]}" run \
  -m "$metadata_file" \
  "$workflow" \
  -i "$inputs_file"

python3 "$stage_script" \
  --metadata "$metadata_file" \
  --sample-list "$sample_list" \
  --output-root "$output_root"

python3 "$repo_root/resources/scripts/py/summarize_timing.py" "$metadata_file" || true

if [[ $cleanup_executions == true ]]; then
  workflow_root=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1])).get("workflowRoot", ""))' "$metadata_file")
  executions_root=$(realpath -m -- "$repo_root/outputs/cromwell-executions")
  if [[ $workflow_root != /* ]]; then
    workflow_root="$repo_root/$workflow_root"
  fi
  workflow_root=$(realpath -m -- "$workflow_root")
  if [[ $workflow_root == "$executions_root"/*/* ]]; then
    rm -rf -- "$workflow_root"
    printf 'Removed execution directory %s\n' "$workflow_root"
  else
    printf 'Skipping cleanup; unexpected workflow root: %s\n' "$workflow_root" >&2
  fi
fi
