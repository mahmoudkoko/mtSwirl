#!/bin/bash
# shellcheck disable=SC2154  # the options are assigned dynamically by `declare`
# Subset a WGS BAM/CRAM to the mt (and optionally NUMT) intervals, drop reads whose mate is
# missing, revert to an unmapped BAM, collect WGS metrics and mean/median coverage, and
# produce a duplicate-marked, coordinate-sorted BAM. Called by MongoSubsetBamToChrMAndRevert
# in MongoTasks.wdl. All outputs are written under ./out (the WDL collects them from there).
#
# Usage: subset_bam_to_chrm_and_revert.sh --sample-name NAME ... (see `required`)
# Optional: --ref-fasta, --mt-interval-list, --nuc-interval-list, --contig-name,
#   --requester-pays-project, --read-name-regex, --coverage-cap, --printreads-extra-args
# --printreads-extra-args is word-split on purpose (it holds extra GATK arguments).

required=(
  sample_name command_mem read_length input_bam input_bai
  skip_restore_hardclips coverage_metrics_script
)
optional=(
  ref_fasta mt_interval_list nuc_interval_list contig_name requester_pays_project
  read_name_regex coverage_cap printreads_extra_args
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

# Not using pipefail: the read-list pipeline below relies on grep returning 1 when no reads match.
set -e

mkdir out

this_sample=out/"${sample_name}"
java_opts="-Xmx${command_mem}m"

printreads_args=()
[[ -n ${ref_fasta:-} ]] && printreads_args+=(-R "$ref_fasta")
[[ -n ${mt_interval_list:-} ]] && printreads_args+=(-L "$mt_interval_list")
[[ -n ${nuc_interval_list:-} ]] && printreads_args+=(-L "$nuc_interval_list")
[[ -n ${contig_name:-} ]] && printreads_args+=(-L "$contig_name")
[[ -n ${requester_pays_project:-} ]] && printreads_args+=(--gcs-project-for-requester-pays "$requester_pays_project")

hardclip_args=()
if [[ $skip_restore_hardclips == true ]]; then
  hardclip_args=(--RESTORE_HARDCLIPS false)
fi

# shellcheck disable=SC2086  # printreads_extra_args is intentionally word-split
gatk --java-options "$java_opts" PrintReads \
  "${printreads_args[@]}" \
  --read-filter MateOnSameContigOrNoMappedMateReadFilter \
  --read-filter MateUnmappedAndUnmappedReadFilter \
  -I "${input_bam}" --read-index "${input_bai}" \
  -O "${this_sample}.bam" ${printreads_extra_args:-}

echo "Now removing mapping..."
set +e
gatk --java-options "$java_opts" ValidateSamFile \
  -INPUT "${this_sample}.bam" \
  -O output.txt \
  -M VERBOSE \
  -IGNORE_WARNINGS true \
  -MAX_OUTPUT 9999999
grep 'ERROR.*Mate not found for paired read' output.txt |
  sed -e 's/ERROR::MATE_NOT_FOUND:Read name //g' |
  sed -e 's/, Mate not found for paired read//g' > read_list.txt
wc -l < read_list.txt | sed 's/^ *//g' > "${this_sample}.ct_failed.txt"
if [[ $(tr -d "\r\n" < read_list.txt | wc -c) -eq 0 ]]; then
  cp "${this_sample}.bam" rescued.bam
else
  gatk --java-options "$java_opts" FilterSamReads \
    -I "${this_sample}.bam" \
    -O rescued.bam \
    -READ_LIST_FILE read_list.txt \
    -FILTER excludeReadList
fi
gatk --java-options "$java_opts" RevertSam \
  -INPUT rescued.bam \
  -OUTPUT_BY_READGROUP false \
  -OUTPUT "${this_sample}.unmap.bam" \
  -VALIDATION_STRINGENCY LENIENT \
  -ATTRIBUTE_TO_CLEAR FT \
  -ATTRIBUTE_TO_CLEAR CO \
  -SORT_ORDER queryname \
  -RESTORE_ORIGINAL_QUALITIES false "${hardclip_args[@]}"

set -e
echo "Now getting WGS metrics on the subsetted bam..."
wgs_args=()
[[ -n ${mt_interval_list:-} ]] && wgs_args+=("INTERVALS=${mt_interval_list}")
[[ -n ${coverage_cap:-} ]] && wgs_args+=("COVERAGE_CAP=${coverage_cap}")
gatk --java-options "$java_opts" CollectWgsMetrics \
  INPUT="${this_sample}.bam" \
  "${wgs_args[@]}" \
  VALIDATION_STRINGENCY=SILENT \
  REFERENCE_SEQUENCE="${ref_fasta:-}" \
  OUTPUT="${this_sample}.wgs_metrics.txt" \
  USE_FAST_ALGORITHM=true \
  READ_LENGTH="${read_length}" \
  INCLUDE_BQ_HISTOGRAM=true \
  THEORETICAL_SENSITIVITY_OUTPUT="${this_sample}.theoretical_sensitivity.txt"

Rscript --vanilla "${coverage_metrics_script}" \
  "${this_sample}.wgs_metrics.txt" \
  "${this_sample}.mean_coverage.txt" \
  "${this_sample}.median_coverage.txt"

echo "Now preprocessing subsetted bam..."
md_args=()
[[ -n ${read_name_regex:-} ]] && md_args+=("READ_NAME_REGEX=${read_name_regex}")
gatk --java-options "$java_opts" MarkDuplicates \
  INPUT="${this_sample}.bam" \
  OUTPUT=md.bam \
  METRICS_FILE="${this_sample}.duplicate.metrics" \
  VALIDATION_STRINGENCY=SILENT \
  "${md_args[@]}" \
  OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
  ASSUME_SORT_ORDER="queryname" \
  CLEAR_DT="false" \
  ADD_PG_TAG_TO_READS=false

gatk --java-options "$java_opts" SortSam \
  INPUT=md.bam \
  OUTPUT="${this_sample}.proc.bam" \
  SORT_ORDER="coordinate" \
  CREATE_INDEX=true \
  MAX_RECORDS_IN_RAM=300000
