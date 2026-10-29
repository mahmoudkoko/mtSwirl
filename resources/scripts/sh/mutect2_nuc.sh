#!/bin/bash
# shellcheck disable=SC2154  # the options are assigned dynamically by `declare`
# Call nuclear variants with Mutect2, filter them (FilterMutectCalls, allele-specific filters,
# optional blacklist mask), keep PASS sites and split multi-allelics. Called by MongoNucM2 in
# MongoTasks.wdl. All outputs are written under ./out (the WDL collects them from there).
#
# Usage: mutect2_nuc.sh --sample-name NAME ... (see `required`)
# Optional: --suffix, --mt-interval-list, --m2-extra-args, --m2-extra-filtering-args,
#   --vaf-filter-threshold, --f-score-beta, --blacklisted-sites
# --m2-extra-args and --m2-extra-filtering-args are word-split on purpose (extra GATK arguments).

required=(
  sample_name command_mem ref_fasta input_bam make_bamout max_reads_per_alignment_start
  max_alt_allele_count nuc_contamination_estimate
)
optional=(
  suffix mt_interval_list m2_extra_args m2_extra_filtering_args vaf_filter_threshold
  f_score_beta blacklisted_sites
)

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

echo "Extra arguments for mutect2: ${m2_extra_args:-}"

mkdir out
this_sample=out/"${sample_name}"
this_basename="${this_sample}${suffix:-}"
bamoutfile="${this_basename}.bamout.bam"
touch "${bamoutfile}"

java_opts="-Xmx${command_mem}m"

m2_args=()
[[ -n ${mt_interval_list:-} ]] && m2_args+=(-L "$mt_interval_list")
[[ -n ${vaf_filter_threshold:-} ]] && m2_args+=(--minimum-allele-fraction "$vaf_filter_threshold")
if [[ $make_bamout == true ]]; then
  m2_args+=(--bam-output "${bamoutfile}")
fi

filter_args=()
[[ -n ${vaf_filter_threshold:-} ]] && filter_args+=(--min-allele-fraction "$vaf_filter_threshold")
[[ -n ${f_score_beta:-} ]] && filter_args+=(--f-score-beta "$f_score_beta")

mask_args=()
[[ -n ${blacklisted_sites:-} ]] && mask_args+=(--mask-name 'blacklisted_site' --mask "$blacklisted_sites")

# shellcheck disable=SC2086  # m2_extra_args is intentionally word-split
gatk --java-options "$java_opts" Mutect2 \
  -R "${ref_fasta}" \
  -I "${input_bam}" \
  -O "${this_basename}.raw.vcf" \
  ${m2_extra_args:-} \
  --annotation StrandBiasBySample \
  --max-reads-per-alignment-start "${max_reads_per_alignment_start}" \
  --max-mnp-distance 0 "${m2_args[@]}"

echo "Now filtering..."
# shellcheck disable=SC2086  # m2_extra_filtering_args is intentionally word-split
gatk --java-options "$java_opts" FilterMutectCalls -V "${this_basename}.raw.vcf" \
  -R "${ref_fasta}" \
  -O filtered.vcf \
  --stats "${this_basename}.raw.vcf.stats" \
  ${m2_extra_filtering_args:-} \
  --max-alt-allele-count "${max_alt_allele_count}" \
  "${filter_args[@]}" \
  --contamination-estimate "${nuc_contamination_estimate}"

gatk --java-options "$java_opts" VariantFiltration -V filtered.vcf \
  -O "${this_basename}.vcf" \
  --apply-allele-specific-filters \
  "${mask_args[@]}"

echo "Now filtering VCF..."
gatk --java-options "$java_opts" SelectVariants \
  -V "${this_basename}.vcf" \
  --exclude-filtered \
  -O "${this_basename}.pass.vcf"

gatk CountVariants -V "${this_basename}.pass.vcf" | tail -n1 > "${this_basename}.passvars.txt"

echo "Now splitting..."
gatk --java-options "$java_opts" LeftAlignAndTrimVariants \
  -R "${ref_fasta}" \
  -V "${this_basename}.pass.vcf" \
  -O "${this_basename}.pass.split.vcf" \
  --split-multi-allelics \
  --dont-trim-alleles \
  --keep-original-ac \
  --create-output-variant-index
