#!/usr/bin/env bash
set -euo pipefail

# Keep all scratch space under the mounted outputs directory so it is
# writable and visible on the host (also under read-only root filesystems).
mkdir -p /opt/mtSwirl/outputs/container-tmp
scratch=$(mktemp -d /opt/mtSwirl/outputs/container-tmp/run.XXXXXX)
trap 'rm -rf -- "$scratch"; rmdir /opt/mtSwirl/outputs/container-tmp 2>/dev/null || true' EXIT

export TMPDIR="$scratch" HOME="$scratch" SPARK_LOCAL_DIRS="$scratch"
export JAVA_TOOL_OPTIONS="-Duser.home=$scratch -Djava.io.tmpdir=$scratch"

# Hadoop/Spark (Hail) fail to start when the UID has no passwd entry, as
# happens with `docker run --user UID:GID`. Fake one for the current user.
if ! getent passwd "$(id -u)" >/dev/null; then
  printf 'mtswirl:x:%s:%s:mtSwirl:%s:/bin/bash\n' "$(id -u)" "$(id -g)" "$scratch" > "$scratch/passwd"
  printf 'mtswirl:x:%s:\n' "$(id -g)" > "$scratch/group"
  export NSS_WRAPPER_PASSWD="$scratch/passwd" NSS_WRAPPER_GROUP="$scratch/group"
  export LD_PRELOAD=libnss_wrapper.so
fi

bash /opt/mtSwirl/run_swirl_pipeline.sh "$@"
