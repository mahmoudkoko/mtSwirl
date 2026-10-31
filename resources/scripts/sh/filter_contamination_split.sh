#!/bin/bash
# shellcheck disable=SC2154  # the options are assigned dynamically by `declare`
# Long-read (Clair3 / ClairS-TO) version. Re-filter the round-1 mt calls: LowAF at the VAF threshold
# (per ALT allele after splitting), blacklist mask, and multi-allelic splitting. The haplocheck
# contamination estimate is recorded for reporting only: Clair has no contamination model.
# Called by MongoM2FilterContaminationSplit in MongoTasks.wdl. All outputs are written under
# ./out (the WDL collects them from there).
#
# Usage: filter_contamination_split.sh --sample-name NAME ... (see `required`)
# Optional: --vaf-filter-threshold, --blacklisted-sites; --m2-extra-filtering-args, --f-score-beta,
#   --raw-vcf-stats and --max-alt-allele-count are accepted but unused (Mutect2 leftovers)

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

set -euo pipefail

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

cp "${raw_vcf}" step.vcf
if [[ -n ${vaf_filter_threshold:-} ]]; then
  bcftools filter -s LowAF -e "MAX(FMT/AF) < ${vaf_filter_threshold}" -Ov -o step_next.vcf step.vcf
  mv step_next.vcf step.vcf
fi
if [[ -n ${blacklisted_sites:-} ]]; then
  bcftools filter -s blacklisted_site -m + --mask-file "${blacklisted_sites}" -Ov -o step_next.vcf step.vcf
  mv step_next.vcf step.vcf
fi
cp step.vcf "${this_basename}.vcf"

echo "Now splitting multi-allelics..."
bcftools norm -m -any -f "${ref_fasta}" -Ov -o split_all.vcf "${this_basename}.vcf"
if [[ -n ${vaf_filter_threshold:-} ]]; then
  # AF floor per ALT allele (a multi-allelic record is only filtered above by its largest ALT).
  bcftools filter -s LowAF -m + -e "FMT/AF < ${vaf_filter_threshold}" -Ov -o "${this_basename}.split.vcf" split_all.vcf
else
  cp split_all.vcf "${this_basename}.split.vcf"
fi
gatk IndexFeatureFile -I "${this_basename}.vcf"
gatk IndexFeatureFile -I "${this_basename}.split.vcf"
