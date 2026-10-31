#!/usr/bin/env python3
"""Print per-stage wall-clock times from a Cromwell metadata JSON."""
import json
import sys
from datetime import datetime


def parse(value):
    return datetime.fromisoformat(value.replace("Z", "+00:00"))


def minutes(call):
    if not call.get("start") or not call.get("end"):
        return None
    return (parse(call["end"]) - parse(call["start"])).total_seconds() / 60


def main():
    with open(sys.argv[1], encoding="utf-8") as handle:
        metadata = json.load(handle)
    print("Stage timing (minutes):")
    for name, calls in metadata.get("calls", {}).items():
        values = [minutes(call) for call in calls]
        values = [value for value in values if value is not None]
        if values:
            label = name.split(".", 1)[-1]
            print(f"  {label:<40} {max(values):7.1f}" + (f"  (x{len(values)} shards, max)" if len(values) > 1 else ""))
    total = minutes(metadata)
    if total is not None:
        print(f"  {'TOTAL':<40} {total:7.1f}")


if __name__ == "__main__":
    main()
