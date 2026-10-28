#!/usr/bin/env bash
set -euo pipefail

usage() {
  printf 'Usage: %s [--validate-only] [--test-run | --sample-manifest FILE] [--configs-json FILE] [--files-json FILE] [--scripts-json FILE]\n' "$0" >&2
  printf 'Use --sample-manifest to specify a CSV, TSV, or TXT sample manifest.\n' >&2
  printf 'Use --test-run to run with the bundled two-sample manifest at inputs/test/samples.tsv.\n' >&2
  printf 'Uses the packaged JSON files configured in inputs/json/pipeline_configs.json by default.\n' >&2
}

validate_only=false
test_run=false
configs_json=""
files_json=""
scripts_json=""
manifest=""
manifest_seen=false

while (($#)); do
  case "$1" in
    --validate-only)
      validate_only=true
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

if [[ $test_run == true && $manifest_seen == true ]]; then
  printf 'Do not combine --test-run with --sample-manifest; --test-run uses the bundled test manifest.\n' >&2
  usage
  exit 2
fi

if [[ $test_run == false && $manifest_seen == false ]]; then
  usage
  exit 2
fi

repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
if [[ $test_run == true ]]; then
  manifest="$repo_root/inputs/test/samples.tsv"
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
if ((${#config_values[@]} != 12)); then
  printf 'Pipeline config did not return all 12 expected values.\n' >&2
  exit 1
fi

hail_path=${config_values[0]}
tools_path=${config_values[1]}
cromwell_jar=${config_values[2]}
womtool_jar=${config_values[3]}
cromwell_config=${config_values[4]}
workflow=${config_values[5]}
prepare_script=${config_values[6]}
stage_script=${config_values[7]}
output_root=${config_values[8]}
files_json=${config_values[9]}
scripts_json=${config_values[10]}
haplocheck_cli=${config_values[11]}

for file in "$cromwell_jar" "$womtool_jar" "$cromwell_config" "$workflow" "$prepare_script" "$stage_script" "$files_json" "$scripts_json" "$haplocheck_cli"; do
  if [[ ! -f $file ]]; then
    printf 'Configured file not found: %s\n' "$file" >&2
    exit 1
  fi
done
configured_bin_dirs() {
  local path=$1
  [[ -d $path/bin ]] && printf '%s\n' "$path/bin"
  [[ -d $path ]] && printf '%s\n' "$path"
}

mapfile -t tools_bin_dirs < <(configured_bin_dirs "$tools_path")
mapfile -t hail_bin_dirs < <(configured_bin_dirs "$hail_path")
path_prefix=""
for bin_dir in "${tools_bin_dirs[@]}" "${hail_bin_dirs[@]}"; do
  path_prefix+="${path_prefix:+:}$bin_dir"
done
if [[ -n $path_prefix ]]; then
  export PATH="$path_prefix:$PATH"
fi

find_executable() {
  local dirs_name=$1
  local name dir
  shift
  local -n dirs_ref=$dirs_name
  for dir in "${dirs_ref[@]}"; do
    for name in "$@"; do
      if [[ -x "$dir/$name" && ! -d "$dir/$name" ]]; then
        printf '%s\n' "$dir/$name"
        return 0
      fi
    done
  done
  return 1
}

java_executable=$(command -v java || true)
if [[ -z $java_executable ]]; then
  printf 'java was not found in tools_path or PATH.\n' >&2
  exit 1
fi
if ! command -v Rscript >/dev/null 2>&1; then
  printf 'Rscript was not found in tools_path or PATH.\n' >&2
  exit 1
fi

python_has_hail() {
  "$1" -c 'import importlib.util, sys; sys.exit(importlib.util.find_spec("hail") is None)' >/dev/null 2>&1
}

hail_python=""
for bin_dir in "${hail_bin_dirs[@]}"; do
  for python_name in python python3; do
    candidate="$bin_dir/$python_name"
    if [[ -x $candidate && ! -d $candidate ]] && python_has_hail "$candidate"; then
      hail_python=$candidate
      break 2
    fi
  done
done
if [[ -z $hail_python ]]; then
  IFS=: read -r -a path_dirs <<< "$PATH"
  for bin_dir in "${path_dirs[@]}"; do
    for python_name in python python3; do
      candidate="$bin_dir/$python_name"
      if [[ -x $candidate && ! -d $candidate ]] && python_has_hail "$candidate"; then
        hail_python=$candidate
        break 2
      fi
    done
  done
fi
if [[ -z $hail_python ]]; then
  printf 'A Python executable with Hail installed was not found in hail_path or PATH.\n' >&2
  exit 1
fi

java_home_for() {
  local java_path resolved_java java_dir
  java_path=$1
  resolved_java=$(readlink -f -- "$java_path")
  java_dir=$(dirname -- "$resolved_java")
  if [[ $(basename -- "$java_dir") == bin ]]; then
    dirname -- "$java_dir"
  else
    printf '%s\n' "$java_dir"
  fi
}

hail_java=$(find_executable hail_bin_dirs java || true)
export JAVA_HOME
JAVA_HOME=$(java_home_for "$java_executable")
if [[ -n $hail_java ]]; then
  HAIL_JAVA_HOME=$(java_home_for "$hail_java")
else
  HAIL_JAVA_HOME=$JAVA_HOME
fi
export HAIL_JAVA_HOME
export HAIL_BIN_DIR
HAIL_BIN_DIR=$(dirname -- "$hail_python")
export HAIL_PYTHON="$hail_python"
export PYSPARK_PYTHON="$hail_python"
export TOOLS_PATH="$tools_path"
export HAIL_PATH="$hail_path"

metadata_dir=$(mktemp -d "${TMPDIR:-/tmp}/mtSwirl-scatter.XXXXXX")
metadata_file="$metadata_dir/metadata.json"
inputs_file="$metadata_dir/inputs.json"
cram_list="$metadata_dir/crams.list"
index_list="$metadata_dir/indexes.list"
sample_list="$metadata_dir/samples.list"
cleanup() {
  rm -f -- "$metadata_file" "$inputs_file" "$cram_list" "$index_list" "$sample_list"
  rmdir -- "$metadata_dir"
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
java -jar "$womtool_jar" validate -i "$inputs_file" -l "$workflow"
if [[ $validate_only == true ]]; then
  printf 'Validation only; Cromwell was not started.\n'
  exit 0
fi

java "-Dconfig.file=$cromwell_config" -jar "$cromwell_jar" run \
  -m "$metadata_file" \
  "$workflow" \
  -i "$inputs_file"

python3 "$stage_script" \
  --metadata "$metadata_file" \
  --sample-list "$sample_list" \
  --output-root "$output_root"
