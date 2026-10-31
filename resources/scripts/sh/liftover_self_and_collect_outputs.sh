#!/bin/bash
# shellcheck disable=SC2154  # the options are assigned dynamically by `declare`
# Lift the per-base self-reference coverage back to GRCh38 (liftOver, then fill in deleted
# sites) and write the per-sample diagnostic statistics table. Called by
# MongoLiftoverSelfAndCollectOutputs in MongoTasks.wdl. All outputs are written under ./out
# (the WDL collects them from there).
#
# Usage: liftover_self_and_collect_outputs.sh --sample-name NAME ... (see `required`)

required=(
  sample_name self_ref_table chain homoplasmic_deletions_coverage liftover_table
  coverage_helpers_script mean_coverage median_coverage major_haplogroup contamination
  nuc_variants_pass n_reads_unpaired_dropped nuc_variants_dropped mtdna_consensus_overlaps
  nuc_consensus_overlaps
)

while [[ $# -gt 0 ]]; do
  key=${1#--}
  key=${key//-/_}
  if [[ $1 != --* || $# -lt 2 || ! " ${required[*]} " == *" $key "* ]]; then
    echo "ERROR: unknown or incomplete option: $1" >&2
    exit 2
  fi
  declare "$key=$2"
  shift 2
done
for key in "${required[@]}"; do
  if [[ -z ${!key+x} ]]; then
    echo "ERROR: missing required option --${key//_/-}" >&2
    exit 2
  fi
done

set -e

mkdir out

this_sample=out/"${sample_name}"

tail -n +2 "${self_ref_table}" | awk '{OFS="\t"} {print $1"\t"$2-1"\t"$2"\t"$4}' > per_base_coverage.bed
liftOver per_base_coverage.bed "${chain}" "${this_sample}.liftedOver.bed" "${this_sample}.liftedOver.unmapped"

Rscript --vanilla "${coverage_helpers_script}" append-lifted-coverage \
  "${chain}" \
  "${this_sample}.liftedOver.unmapped" \
  "${this_sample}.liftedOver.bed" \
  "${homoplasmic_deletions_coverage}" \
  "${this_sample}.appended.liftedOver.tsv"

Rscript --vanilla "${coverage_helpers_script}" write-diagnostic-statistics \
  "${liftover_table}" \
  "${this_sample}_mtanalysis_diagnostic_statistics.tsv" \
  "${mean_coverage}" \
  "${median_coverage}" \
  "${major_haplogroup}" \
  "${contamination}" \
  "${nuc_variants_pass}" \
  "${n_reads_unpaired_dropped}" \
  "${nuc_variants_dropped}" \
  "${mtdna_consensus_overlaps}" \
  "${nuc_consensus_overlaps}"
