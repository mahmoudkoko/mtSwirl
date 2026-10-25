#! /bin/bash

gcloud storage cp \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.dict" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.amb" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.ann" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.bwt" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.fai" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.pac" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.sa" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.dict" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.amb" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.ann" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.bwt" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.fai" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.pac" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.sa" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/README" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/ShiftBack.chain" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/blacklist_sites.hg38.chrM.bed" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/blacklist_sites.hg38.chrM.bed.idx" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/blacklist_sites.hg38.chrM.shifted_by_8000_bases.bed.idx.old" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/blacklist_sites.hg38.chrM.shifted_by_8000_bases.bed.old" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/blacklist_sites.hg38.chrM.shifted_by_8000_bases.fixed.bed" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/blacklist_sites.hg38.chrM.shifted_by_8000_bases.fixed.bed.idx" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/chrM.hg38.interval_list" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/chrMWithFinalNuMTs.hg38.interval_list" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/control_region_shifted.chrM.interval_list" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/non_control_region.chrM.interval_list" \
  .
