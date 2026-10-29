#!/bin/bash
# shellcheck disable=SC2154  # the options are assigned dynamically by `declare`
# Re-filter the mt Mutect2 calls using the larger of the nuclear and haplocheck contamination
# estimates, apply allele-specific filters and the blacklist mask, and split multi-allelics.
# Called by MongoM2FilterContaminationSplit in MongoTasks.wdl. All outputs are written under
# ./out (the WDL collects them from there).
#
# Usage: filter_contamination_split.sh --sample-name NAME ... (see `required`)
# Optional: --m2-extra-filtering-args (word-split on purpose), --vaf-filter-threshold,
#   --f-score-beta, --blacklisted-sites

required=(
  sample_name suffix command_mem ref_fasta raw_vcf raw_vcf_stats has_contamination
  contamination_major contamination_minor nuc_contamination_estimate max_alt_allele_count
)
optional=(m2_extra_filtering_args vaf_filter_threshold f_score_beta blacklisted_sites)

while [[ $# -gt 0 ]]; do
  key=${1#--}
  key=${key//-/_}
  if [[ $1 != --* || $# -lt 2 || ! " ${required[*]} ${optional[*]} " == *" $key "* ]]; then
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
this_basename="${this_sample}${suffix}"
touch "${this_basename}.bamout.bam"

if [[ "${has_contamination}" == 'YES' ]]; then
  if (( $(echo "${contamination_major} == 0.0" | bc -l) )); then
    this_hc_contamination="${contamination_minor}"
  else
    this_hc_contamination=$(bc <<< "1-${contamination_major}")
  fi
else
  this_hc_contamination=0.0
fi

echo "${this_hc_contamination}" > "${this_basename}.hc_contam.txt"

if (( $(echo "${nuc_contamination_estimate} > ${this_hc_contamination}" | bc -l) )); then
  this_max_contamination="${nuc_contamination_estimate}"
else
  this_max_contamination="${this_hc_contamination}"
fi

echo "Nuclear estimate: ${nuc_contamination_estimate}; HC: ${this_hc_contamination}; Max: ${this_max_contamination}"

java_opts="-Xmx${command_mem}m"

filter_args=()
[[ -n ${vaf_filter_threshold:-} ]] && filter_args+=(--min-allele-fraction "$vaf_filter_threshold")
[[ -n ${f_score_beta:-} ]] && filter_args+=(--f-score-beta "$f_score_beta")

mask_args=()
[[ -n ${blacklisted_sites:-} ]] && mask_args+=(--mask-name 'blacklisted_site' --mask "$blacklisted_sites")

# to do: check why this is fixed not dynamic
# shellcheck disable=SC2086  # m2_extra_filtering_args is intentionally word-split
gatk --java-options "$java_opts" FilterMutectCalls \
  -V "${raw_vcf}" \
  -R "${ref_fasta}" \
  -O filtered.vcf \
  --stats "${raw_vcf_stats}" \
  ${m2_extra_filtering_args:-} \
  --max-alt-allele-count "${max_alt_allele_count}" \
  --mitochondria-mode \
  "${filter_args[@]}" \
  --contamination-estimate "${this_max_contamination}"
# The static blacklist BED and its GATK index are supplied as workflow inputs.
# Do not regenerate the index here: Cromwell may hard-link localized inputs.

gatk --java-options "$java_opts" VariantFiltration \
  -V filtered.vcf \
  -O "${this_basename}.vcf" \
  --apply-allele-specific-filters \
  "${mask_args[@]}"

echo "Now splitting multi-allelics..."
gatk --java-options "$java_opts" LeftAlignAndTrimVariants \
  -R "${ref_fasta}" \
  -V "${this_basename}.vcf" \
  -O "${this_basename}.split.vcf" \
  --split-multi-allelics \
  --dont-trim-alleles \
  --keep-original-ac \
  --create-output-variant-index
