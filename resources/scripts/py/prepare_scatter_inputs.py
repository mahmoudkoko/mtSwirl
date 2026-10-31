#!/usr/bin/env python3

import argparse
import csv
import json
import re
import sys
from pathlib import Path


ALIASES = {
    "sample": {"sample", "sampleid", "samplename", "id"},
    "bam": {
        "bam",
        "bamfile",
        "bampath",
        "cram",
        "cramfile",
        "crampath",
        "wgsalignedinputbamorcram",
    },
    "index": {
        "index",
        "indexfile",
        "indexpath",
        "idx",
        "idxfile",
        "idxpath",
        "bai",
        "crai",
        "wgsalignedinputbamorcramindex",
    },
}
WRAPPER_PREFIX = "SwirlPipeline."
PIPELINE_PREFIX = "MitochondriaPipeline."
HELPER_SCRIPT_KEYS = {
    "SwirlPipeline.prepare_script",
    "SwirlPipeline.stage_script",
}


def parse_args():
    parser = argparse.ArgumentParser(description="Prepare scatter workflow input JSON.")
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--repo-root", required=True, type=Path)
    parser.add_argument("--files-json", required=True, type=Path)
    parser.add_argument("--scripts-json", required=True, type=Path)
    parser.add_argument("--configs-json", required=True, type=Path)
    parser.add_argument("--haplocheck-cli", required=True, type=Path)
    parser.add_argument("--cram-list", required=True, type=Path)
    parser.add_argument("--index-list", required=True, type=Path)
    parser.add_argument("--sample-list", required=True, type=Path)
    parser.add_argument("--inputs-output", required=True, type=Path)
    return parser.parse_args()


def normalize_header(value):
    return re.sub(r"[^a-z0-9]", "", value.strip().lower())


def resolve_manifest_path(value, manifest, repo_root, row_number, column):
    path = Path(value.strip()).expanduser()
    if not path.is_absolute():
        candidates = (manifest.parent / path, repo_root / path)
        path = next((candidate for candidate in candidates if candidate.is_file()), None)
        if path is None:
            choices = " or ".join(str(candidate.resolve()) for candidate in candidates)
            raise ValueError(f"Row {row_number}: {column} file does not exist: {value!r} (checked {choices})")
    if not path.is_file():
        raise ValueError(f"Row {row_number}: {column} is not a file: {path}")
    return path.resolve()


def read_manifest(manifest, repo_root):
    with manifest.open(newline="", encoding="utf-8-sig") as source:
        sample = source.read(8192)
        source.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        except csv.Error as error:
            raise ValueError(f"Could not determine whether the manifest is CSV, TSV, or delimited text: {error}") from error

        reader = csv.reader(source, dialect)
        try:
            headers = next(reader)
        except StopIteration as error:
            raise ValueError("The sample manifest is empty.") from error

        column_indexes = {}
        for kind, choices in ALIASES.items():
            matches = [i for i, header in enumerate(headers) if normalize_header(header) in choices]
            if len(matches) != 1:
                expected = ", ".join(sorted(choices))
                raise ValueError(f"Manifest must have exactly one {kind} column (accepted headers: {expected}).")
            column_indexes[kind] = matches[0]

        samples, crams, indexes = [], [], []
        seen_samples = set()
        for row_number, row in enumerate(reader, start=2):
            if not row or all(not value.strip() for value in row):
                continue
            if len(row) != len(headers):
                raise ValueError(f"Row {row_number}: expected {len(headers)} columns, found {len(row)}.")
            values = {kind: row[index].strip() for kind, index in column_indexes.items()}
            sample_id = values["sample"]
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", sample_id):
                raise ValueError(f"Row {row_number}: invalid sample ID {sample_id!r}.")
            if sample_id in seen_samples:
                raise ValueError(f"Row {row_number}: duplicate sample ID {sample_id!r}.")
            if any(not values[kind] for kind in ("bam", "index")):
                raise ValueError(f"Row {row_number}: BAM/CRAM and index paths must not be empty.")
            seen_samples.add(sample_id)
            samples.append(sample_id)
            crams.append(resolve_manifest_path(values["bam"], manifest, repo_root, row_number, "BAM/CRAM"))
            indexes.append(resolve_manifest_path(values["index"], manifest, repo_root, row_number, "index"))

    if not samples:
        raise ValueError("The sample manifest contains no data rows.")
    return samples, crams, indexes


def wrapper_key(key):
    if key.startswith(WRAPPER_PREFIX):
        return key
    if key.startswith(PIPELINE_PREFIX):
        return WRAPPER_PREFIX + key[len(PIPELINE_PREFIX) :]
    raise ValueError(
        f"Input JSON key {key!r} must begin with {PIPELINE_PREFIX!r} "
        f"or {WRAPPER_PREFIX!r}."
    )


def read_json_inputs(path, repo_root, resolve_files):
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(f"Invalid JSON in {path}: {error}") from error
    if not isinstance(data, dict):
        raise ValueError(f"Input JSON must contain a JSON object: {path}")

    normalized = {}
    for key, value in data.items():
        if not isinstance(key, str):
            raise ValueError(f"Input JSON keys must be strings: {path}")
        if key in HELPER_SCRIPT_KEYS:
            continue
        output_key = wrapper_key(key)
        if resolve_files:
            if not isinstance(value, str) or not value:
                raise ValueError(f"File input {key!r} in {path} must have a non-empty path string.")
            file_path = Path(value).expanduser()
            if not file_path.is_absolute():
                file_path = repo_root / file_path
            if not file_path.is_file():
                raise ValueError(f"File input {key!r} does not exist: {file_path}")
            value = str(file_path.resolve())
        merge_inputs(normalized, {output_key: value}, path)
    return normalized


def merge_inputs(target, incoming, source):
    for key, value in incoming.items():
        if key in target and target[key] != value:
            raise ValueError(f"Conflicting values for {key!r} while loading {source}.")
        target[key] = value


def read_json_options(data, source):
    normalized = {}
    for key, value in data.items():
        if not isinstance(key, str):
            raise ValueError(f"Workflow option keys must be strings: {source}")
        output_key = wrapper_key(key)
        merge_inputs(normalized, {output_key: value}, source)
    return normalized


def main():
    args = parse_args()
    repo_root = args.repo_root.resolve()
    manifest = args.manifest.resolve()
    samples, crams, indexes = read_manifest(manifest, repo_root)

    inputs = {}
    for path in (args.files_json, args.scripts_json):
        merge_inputs(inputs, read_json_inputs(path, repo_root, resolve_files=True), path)
    try:
        configs = json.loads(args.configs_json.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(f"Invalid JSON in {args.configs_json}: {error}") from error
    workflow = configs.get("workflow") if isinstance(configs, dict) else None
    if not isinstance(workflow, dict) or not isinstance(workflow.get("options"), dict):
        raise ValueError("Pipeline config JSON must contain a 'workflow.options' object.")
    merge_inputs(
        inputs,
        read_json_options(workflow["options"], args.configs_json),
        args.configs_json,
    )
    merge_inputs(
        inputs,
        {"SwirlPipeline.haplocheck_cli": str(args.haplocheck_cli.resolve())},
        args.configs_json,
    )

    generated_values = {
        "wgs_aligned_input_bam_or_cram_list": str(args.cram_list.resolve()),
        "wgs_aligned_input_bam_or_cram_index_list": str(args.index_list.resolve()),
        "sample_name_list": str(args.sample_list.resolve()),
    }
    for name in generated_values:
        key = WRAPPER_PREFIX + name
        if key in inputs:
            raise ValueError(f"Input JSONs must not define manifest-generated input {key!r}.")
    inputs.update({WRAPPER_PREFIX + name: value for name, value in generated_values.items()})

    args.cram_list.write_text("".join(f"{path}\n" for path in crams), encoding="utf-8")
    args.index_list.write_text("".join(f"{path}\n" for path in indexes), encoding="utf-8")
    args.sample_list.write_text("".join(f"{sample}\n" for sample in samples), encoding="utf-8")
    args.inputs_output.write_text(json.dumps(inputs, indent=2) + "\n", encoding="utf-8")
    print(f"Prepared {len(samples)} samples from {manifest}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        print(error, file=sys.stderr)
        raise SystemExit(1) from error
