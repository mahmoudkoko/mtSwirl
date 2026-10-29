#!/bin/bash
# shellcheck disable=SC2154  # the options are assigned dynamically by `declare`
# Call nuclear variants with HaplotypeCaller, apply hard filters to SNPs and indels, keep PASS
# sites and split multi-allelics. Called by MongoHC in MongoTasks.wdl. All outputs are written
# under ./out (the WDL collects them from there).
#
# Usage: haplotype_caller_nuc.sh --sample-name NAME ... (see `required`)
# Optional: --suffix, --nuc-interval-list, --hc-extra-args, --force-call-vcf, --contamination
# --hc-extra-args is word-split on purpose (it holds extra GATK arguments).

required=(
  sample_name command_mem ref_fasta input_bam make_bamout
  max_reads_per_alignment_start haplotype_caller_nucdna_dp_lower_bound
)
optional=(suffix nuc_interval_list hc_extra_args force_call_vcf contamination)

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
this_basename="${this_sample}${suffix:-}"
bamoutfile="${this_basename}.bamout.bam"
# Cromwell doesn't like optional task outputs, so we have to touch this file
touch "${bamoutfile}"

java_opts="-Xmx${command_mem}m"

hc_args=()
[[ -n ${nuc_interval_list:-} ]] && hc_args+=(-L "$nuc_interval_list")
[[ -n ${force_call_vcf:-} ]] && hc_args+=(--genotype-filtered-alleles --alleles "$force_call_vcf")
bamout_args=()
if [[ $make_bamout == true ]]; then
  bamout_args=(--bam-output "${bamoutfile}")
fi

dp_filter_args=(
  --genotype-filter-expression "DP < ${haplotype_caller_nucdna_dp_lower_bound}"
  --genotype-filter-name "genoDP${haplotype_caller_nucdna_dp_lower_bound}"
)

# shellcheck disable=SC2086  # hc_extra_args is intentionally word-split
gatk --java-options "$java_opts" HaplotypeCaller \
  -R "${ref_fasta}" \
  -I "${input_bam}" \
  "${hc_args[@]}" \
  -O "${this_basename}.raw.vcf" \
  ${hc_extra_args:-} \
  -contamination "${contamination:-0}" \
  --max-reads-per-alignment-start "${max_reads_per_alignment_start}" \
  --max-mnp-distance 0 \
  --annotation StrandBiasBySample \
  -G StandardAnnotation -G StandardHCAnnotation \
  -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
  "${bamout_args[@]}"

echo "Now applying hard filters..."
gatk --java-options "$java_opts" SelectVariants -V "${this_basename}.raw.vcf" -select-type SNP -O snps.vcf
gatk --java-options "$java_opts" VariantFiltration -V snps.vcf \
  -R "${ref_fasta}" \
  -O snps_filtered.vcf \
  -filter "QD < 2.0" --filter-name "QD2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "SOR > 3.0" --filter-name "SOR3" \
  -filter "FS > 60.0" --filter-name "FS60" \
  -filter "MQ < 40.0" --filter-name "MQ40" \
  -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  --genotype-filter-expression "isHet == 1" --genotype-filter-name "isHetFilt" \
  --genotype-filter-expression "isHomRef == 1" --genotype-filter-name "isHomRefFilt" \
  "${dp_filter_args[@]}"

gatk --java-options "$java_opts" SelectVariants -V "${this_basename}.raw.vcf" -select-type INDEL -O indels.vcf
gatk --java-options "$java_opts" VariantFiltration -V indels.vcf \
  -R "${ref_fasta}" \
  -O indels_filtered.vcf \
  -filter "QD < 2.0" --filter-name "QD2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "FS > 200.0" --filter-name "FS200" \
  -filter "SOR > 10.0" --filter-name "SOR10" \
  -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
  --genotype-filter-expression "isHet == 1" --genotype-filter-name "isHetFilt" \
  --genotype-filter-expression "isHomRef == 1" --genotype-filter-name "isHomRefFilt" \
  "${dp_filter_args[@]}"

gatk --java-options "$java_opts" MergeVcfs -I snps_filtered.vcf -I indels_filtered.vcf -O "${this_basename}.vcf"

echo "Now filtering VCF..."
gatk --java-options "$java_opts" SelectVariants \
  -V "${this_basename}.vcf" \
  --exclude-filtered \
  --set-filtered-gt-to-nocall \
  --exclude-non-variants \
  -O "${this_basename}.pass.vcf"

gatk --java-options "$java_opts" CountVariants -V "${this_basename}.pass.vcf" | tail -n1 > "${this_basename}.passvars.txt"

echo "Now splitting multi-allelics..."
gatk --java-options "$java_opts" LeftAlignAndTrimVariants \
  -R "${ref_fasta}" \
  -V "${this_basename}.pass.vcf" \
  -O "${this_basename}.pass.split.vcf" \
  --split-multi-allelics \
  --dont-trim-alleles \
  --keep-original-ac \
  --create-output-variant-index
