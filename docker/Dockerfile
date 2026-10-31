# Self-contained runtime for the mtSwirl pipeline.
# Cromwell, Womtool and Hail are installed in the single tools environment.
# References are mounted, not baked in.
FROM mambaorg/micromamba:2.0-debian12-slim AS build

USER root
COPY docker/tools.yml /tmp/
RUN micromamba create -y -p /opt/envs/mtSwirl-tools -f /tmp/tools.yml \
 && /opt/envs/mtSwirl-tools/bin/python -m pip install --no-cache-dir "hail==0.2.139" \
 && micromamba clean -a -y \
 && find /opt/envs -name '__pycache__' -type d -prune -exec rm -rf {} +

FROM debian:trixie-slim
RUN apt-get update \
 && apt-get install -y --no-install-recommends bash bc gawk sed grep procps ca-certificates libgomp1 libnss-wrapper \
 && rm -rf /var/lib/apt/lists/*
COPY --from=build /opt/envs /opt/envs
ENV TOOLS_PATH=/opt/envs/mtSwirl-tools \
    PATH=/opt/envs/mtSwirl-tools/bin:/usr/local/bin:/usr/bin:/bin

# Build from the repo root: docker build -f docker/Dockerfile -t mtswirl:mitocamb .
# Mount the reference at /opt/mtSwirl/resources/grch38 and a writable
# directory at /opt/mtSwirl/outputs.
COPY . /opt/mtSwirl
WORKDIR /opt/mtSwirl
ENTRYPOINT ["bash", "/opt/mtSwirl/docker/entrypoint.sh"]
