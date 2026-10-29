#!/bin/bash
# shellcheck disable=SC2154  # the options are assigned dynamically by `declare`
# Realign reads from an unmapped BAM to the personalized mt reference and to its shifted
# version (BWA, MergeBamAlignment, MarkDuplicates, SortSam, ReorderSam), then collect WGS
# metrics and mean/median mt coverage. Called by MongoAlignToMtRegShiftedAndMetrics in
# MongoTasks.wdl. All outputs are written under ./out (the WDL collects them from there).
#
# Usage: align_to_mt_ref_shifted.sh --sample-base-name NAME --suffix S ... (see `required`)
# Optional: --read-name-regex, --coverage-cap
# --bwa-commandline may contain the literal text $bash_ref_fasta, which is replaced by the
# reference being aligned to.

required=(
  sample_base_name suffix command_mem read_length bwa_commandline
  input_bam mt_cat mt_shifted_cat mt mt_dict mt_shifted_dict mt_interval_list
  coverage_metrics_script
)
optional=(read_name_regex coverage_cap)

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

BWAVERSION=$(bwa 2>&1 | grep -e '^Version' | sed 's/Version: //')
export BWAVERSION

set -o pipefail
set -e

mkdir out

this_sample=out/"${sample_base_name}${suffix}"
this_output_bam_basename=out/"$(basename "${input_bam}" .bam).remap${suffix}"

read_name_regex_args=()
if [[ -n ${read_name_regex:-} ]]; then
  read_name_regex_args=(--READ_NAME_REGEX "$read_name_regex")
fi
coverage_cap_args=()
if [[ -n ${coverage_cap:-} ]]; then
  coverage_cap_args=(--COVERAGE_CAP "$coverage_cap")
fi

# Align reads.fastq to one reference and write <base><tag>.bam (reordered to the mt dictionary)
# plus the <base><tag>_pre_mt_filt.bam, BWA stderr log and duplicate metrics.
# tag is empty for the unshifted reference and ".shifted" for the shifted one.
align_and_clean() {
  local tag=$1 cat_fasta=$2 dict=$3
  local bwa_cmd=${bwa_commandline//\$bash_ref_fasta/$cat_fasta}

  bwa index "$cat_fasta"

  # shellcheck disable=SC2086  # the command line is intentionally word-split
  $bwa_cmd reads.fastq - \
    > "aligned${tag}.sam" \
    2> >(tee "${this_output_bam_basename}${tag}.bwa.stderr.log" >&2)

  gatk --java-options "-Xmx${command_mem}m" MergeBamAlignment \
    --VALIDATION_STRINGENCY SILENT \
    --EXPECTED_ORIENTATIONS FR \
    --ATTRIBUTES_TO_RETAIN X0 \
    --ATTRIBUTES_TO_REMOVE NM \
    --ATTRIBUTES_TO_REMOVE MD \
    --ALIGNED_BAM "aligned${tag}.sam" \
    --UNMAPPED_BAM "${input_bam}" \
    -O "mba${tag}.bam" \
    -R "$cat_fasta" \
    --PAIRED_RUN true \
    --SORT_ORDER unsorted \
    --IS_BISULFITE_SEQUENCE false \
    --ALIGNED_READS_ONLY false \
    --CLIP_ADAPTERS false \
    --MAX_RECORDS_IN_RAM 2000000 \
    --ADD_MATE_CIGAR true \
    --MAX_INSERTIONS_OR_DELETIONS -1 \
    --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
    --PROGRAM_RECORD_ID bwamem \
    --PROGRAM_GROUP_VERSION "$BWAVERSION" \
    --PROGRAM_GROUP_COMMAND_LINE "$bwa_cmd" \
    --PROGRAM_GROUP_NAME bwamem \
    --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
    --ALIGNER_PROPER_PAIR_FLAGS true \
    --UNMAP_CONTAMINANT_READS true \
    --ADD_PG_TAG_TO_READS false

  gatk --java-options "-Xmx${command_mem}m" MarkDuplicates \
    -I "mba${tag}.bam" \
    -O "md${tag}.bam" \
    -M "${this_output_bam_basename}${tag}.metrics" \
    --VALIDATION_STRINGENCY SILENT \
    "${read_name_regex_args[@]}" \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
    --ASSUME_SORT_ORDER queryname \
    --CLEAR_DT false \
    --ADD_PG_TAG_TO_READS false

  gatk --java-options "-Xmx${command_mem}m" SortSam \
    -I "md${tag}.bam" \
    -O "${this_output_bam_basename}${tag}_pre_mt_filt.bam" \
    --SORT_ORDER coordinate \
    --CREATE_INDEX true \
    --MAX_RECORDS_IN_RAM 300000

  # now we have to subset to mito and update sequence dictionary
  gatk --java-options "-Xmx${command_mem}m" ReorderSam \
    -I "${this_output_bam_basename}${tag}_pre_mt_filt.bam" \
    -O "${this_output_bam_basename}${tag}.bam" \
    --SEQUENCE_DICTIONARY "$dict" \
    --ALLOW_INCOMPLETE_DICT_CONCORDANCE true \
    --CREATE_INDEX true
}

gatk --java-options "-Xmx${command_mem}m" SamToFastq \
  -I "${input_bam}" \
  -F reads.fastq \
  --INTERLEAVE true \
  --INCLUDE_NON_PF_READS true

align_and_clean "" "$mt_cat" "$mt_dict"

echo "Now starting on shifted..."
align_and_clean ".shifted" "$mt_shifted_cat" "$mt_shifted_dict"

echo "Now collecting wgs metrics..."
gatk --java-options "-Xmx${command_mem}m" CollectWgsMetrics \
  -I "${this_output_bam_basename}.bam" \
  --INTERVALS "${mt_interval_list}" \
  --VALIDATION_STRINGENCY SILENT \
  -R "${mt}" \
  -O "${this_sample}_r2_wgs_metrics.txt" \
  --USE_FAST_ALGORITHM true \
  --READ_LENGTH "${read_length}" \
  "${coverage_cap_args[@]}" \
  --INCLUDE_BQ_HISTOGRAM true \
  --THEORETICAL_SENSITIVITY_OUTPUT "${this_sample}_r2_wgs_theoretical_sensitivity.txt"

Rscript --vanilla "${coverage_metrics_script}" \
  "${this_sample}_r2_wgs_metrics.txt" \
  "${this_sample}_r2_mean_coverage.txt" \
  "${this_sample}_r2_median_coverage.txt"
