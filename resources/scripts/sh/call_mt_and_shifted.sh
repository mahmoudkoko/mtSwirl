#!/bin/bash
# shellcheck disable=SC2154  # the options are assigned dynamically by `declare`
# Call mt variants with Mutect2 (mitochondria mode, force-calling a given VCF) on both the
# regular and the shifted self-reference. Called by MongoCallMtAndShifted in MongoTasks.wdl.
# All outputs are written under ./out (the WDL collects them from there).
#
# Usage: call_mt_and_shifted.sh --sample-base-name NAME --suffix S ... (see `required`)
# Optional: --m2-extra-args, --shifted-m2-extra-args (word-split on purpose: extra GATK arguments)

required=(
  sample_base_name suffix command_mem mem mean_coverage make_bamout max_reads_per_alignment_start
  mt_self input_bam mt_interval_list force_call_vcf
  shifted_mt_self shifted_input_bam shifted_mt_interval_list shifted_force_call_vcf
)
optional=(m2_extra_args shifted_m2_extra_args)

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

# Heap is sized by mem; very deep mitochondrial coverage may need more (assumes ~4000x of coverage per GB).
if [ "${mean_coverage}" -gt "$(( mem * 4000 ))" ]; then
  echo "WARNING: mean mtDNA coverage ${mean_coverage} exceeds the level mem=${mem} GB is expected to handle; raise mem_mutect_mt if Mutect2 runs out of memory."
fi
echo "Extra arguments for mutect2: ${m2_extra_args:-}"

mkdir out

this_sample=out/"${sample_base_name}${suffix}"

touch "${this_sample}.bamout.bam"
touch "${this_sample}.shifted.bamout.bam"

# call_mutect FASTA BAM INTERVALS FORCE_VCF LOCAL_VCF_NAME OUT_PREFIX BAMOUT_FILE EXTRA_ARGS
call_mutect() {
  local fasta=$1 bam=$2 intervals=$3 force_vcf=$4 local_vcf=$5 out_prefix=$6 bamout=$7 extra_args=$8
  local bamout_args=()
  if [[ $make_bamout == true ]]; then
    bamout_args=(--bam-output "$bamout")
  fi

  echo "Obtaining force calls for specified VCF: ${force_vcf}"
  # Cromwell inputs may be read-only symlinks; copy locally before GATK writes the index.
  cp "${force_vcf}" "${local_vcf}"
  gatk IndexFeatureFile -I "${local_vcf}"

  # shellcheck disable=SC2086  # extra_args is intentionally word-split
  gatk --java-options "-Xmx${command_mem}m" Mutect2 \
    -R "${fasta}" \
    -I "${bam}" \
    -L "${intervals}" \
    -O "${out_prefix}.raw.vcf" \
    --genotype-filtered-alleles \
    --alleles "${local_vcf}" \
    ${extra_args} \
    --annotation StrandBiasBySample \
    --mitochondria-mode \
    --read-filter MateOnSameContigOrNoMappedMateReadFilter \
    --read-filter MateUnmappedAndUnmappedReadFilter \
    --max-reads-per-alignment-start "${max_reads_per_alignment_start}" \
    --max-mnp-distance 0 "${bamout_args[@]}"
}

call_mutect "${mt_self}" "${input_bam}" "${mt_interval_list}" "${force_call_vcf}" \
  force_call.vcf.bgz "${this_sample}" "${this_sample}.bamout.bam" "${m2_extra_args:-}"

call_mutect "${shifted_mt_self}" "${shifted_input_bam}" "${shifted_mt_interval_list}" "${shifted_force_call_vcf}" \
  shifted_force_call.vcf.bgz "${this_sample}.shifted" "${this_sample}.shifted.bamout.bam" "${shifted_m2_extra_args:-}"
