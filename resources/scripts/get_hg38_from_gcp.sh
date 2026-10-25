#! /bin/bash

gcloud storage cp \
  "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict" \
  "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta" \
  "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai" \
  .
