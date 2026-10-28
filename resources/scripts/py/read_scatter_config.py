#!/usr/bin/env python3

import argparse
import json
import os
import sys
from pathlib import Path


RUNNER_KEYS = ("cromwell_jar", "womtool_jar", "cromwell_config")
WORKFLOW_KEYS = ("wdl", "files", "scripts", "output_root")
SWIRL_SCRIPTS_PREFIX = "SwirlPipeline."
HELPER_SCRIPTS = ("prepare_script", "stage_script")


def parse_args():
    parser = argparse.ArgumentParser(description="Read and resolve scatter runner configuration.")
    parser.add_argument("--config-json", required=True, type=Path)
    parser.add_argument("--repo-root", required=True, type=Path)
    parser.add_argument("--files-json", type=Path)
    parser.add_argument("--scripts-json", type=Path)
    return parser.parse_args()


def expand_path(value, repo_root):
    expanded = Path(os.path.expandvars(os.path.expanduser(value)))
    if not expanded.is_absolute():
        expanded = repo_root / expanded
    return expanded.resolve()


def require_object(parent, key, source):
    value = parent.get(key)
    if not isinstance(value, dict):
        raise ValueError(f"{source} must contain a '{key}' object.")
    return value


def require_path(mapping, key, section):
    value = mapping.get(key)
    if not isinstance(value, str) or not value:
        raise ValueError(f"Pipeline config '{section}' must define a non-empty path for {key!r}.")
    return value


def optional_path(mapping, key, environment_variable, repo_root):
    value = os.environ.get(environment_variable, mapping.get(key, ""))
    if value is None:
        return ""
    if not isinstance(value, str):
        raise ValueError(f"Pipeline config path {key!r} must be a string.")
    if not value:
        return ""
    return str(expand_path(value, repo_root))


def read_scripts(path, repo_root):
    try:
        scripts = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(f"Invalid JSON in {path}: {error}") from error
    if not isinstance(scripts, dict):
        raise ValueError(f"Scripts JSON must contain an object: {path}")

    resolved = {}
    for name in HELPER_SCRIPTS:
        key = SWIRL_SCRIPTS_PREFIX + name
        value = scripts.get(key)
        if not isinstance(value, str) or not value:
            raise ValueError(f"Scripts JSON must define a non-empty path for {key!r}.")
        script_path = expand_path(value, repo_root)
        if not script_path.is_file():
            raise ValueError(f"Configured helper script does not exist: {script_path}")
        resolved[name] = str(script_path)
    return resolved


def main():
    args = parse_args()
    repo_root = args.repo_root.resolve()
    try:
        config = json.loads(args.config_json.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(f"Invalid JSON in {args.config_json}: {error}") from error
    if not isinstance(config, dict):
        raise ValueError("Pipeline config JSON must contain an object.")
    runner = require_object(config, "runner", "Pipeline config")
    workflow = require_object(config, "workflow", "Pipeline config")
    dependencies = require_object(workflow, "dependencies", "Pipeline config workflow")
    for key in RUNNER_KEYS:
        require_path(runner, key, "runner")
    for key in WORKFLOW_KEYS:
        require_path(workflow, key, "workflow")
    require_path(dependencies, "haplocheck_cli", "workflow.dependencies")
    tools_path = optional_path(dependencies, "tools_path", "TOOLS_PATH", repo_root)
    hail_path = optional_path(dependencies, "hail_path", "HAIL_PATH", repo_root)
    if not isinstance(workflow.get("options"), dict):
        raise ValueError("Pipeline config 'workflow' must contain an 'options' object.")

    files_json = args.files_json or expand_path(workflow["files"], repo_root)
    scripts_json = args.scripts_json or expand_path(workflow["scripts"], repo_root)
    for path in (files_json, scripts_json):
        if not path.is_file():
            raise ValueError(f"Configured input JSON file does not exist: {path}")

    resolved = {}
    for key in RUNNER_KEYS:
        path = expand_path(runner[key], repo_root)
        if not path.is_file():
            raise ValueError(f"Configured runner file does not exist: {path}")
        resolved[key] = str(path)
    for key in WORKFLOW_KEYS:
        path = expand_path(workflow[key], repo_root)
        if key != "output_root" and not path.is_file():
            raise ValueError(f"Configured workflow file does not exist: {path}")
        resolved[key] = str(path)
    resolved.update(read_scripts(scripts_json, repo_root))
    haplocheck_cli = expand_path(dependencies["haplocheck_cli"], repo_root)
    if not haplocheck_cli.is_file():
        raise ValueError(f"Configured HaploCheck JAR does not exist: {haplocheck_cli}")

    values = [
        hail_path,
        tools_path,
        resolved["cromwell_jar"],
        resolved["womtool_jar"],
        resolved["cromwell_config"],
        resolved["wdl"],
        resolved["prepare_script"],
        resolved["stage_script"],
        resolved["output_root"],
        str(files_json.resolve()),
        str(scripts_json.resolve()),
        str(haplocheck_cli),
    ]
    if any("\n" in value or "\t" in value for value in values):
        raise ValueError("Configured paths may not contain tabs or newlines.")
    print("\n".join(values))


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        print(error, file=sys.stderr)
        raise SystemExit(1) from error
