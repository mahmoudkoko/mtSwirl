#!/bin/bash
# shellcheck disable=SC2154  # the options are assigned dynamically by `declare`
# Long-read (ONT) version. Realign reads from an unmapped BAM (modified-base tags kept) to the
# personalized mt reference and to its shifted version with minimap2, filter reads by aligned length
# and aligned fraction (MitoHPC2 filter.lr.sh), restrict to the mt contigs, then collect WGS metrics
# and mean/median mt coverage. Called by MongoAlignToMtRegShiftedAndMetrics in MongoTasks.wdl.
# All outputs are written under ./out (the WDL collects them from there).
#
# Usage: align_to_mt_ref_shifted.sh --sample-base-name NAME --suffix S ... (see `required`)
# Optional: --coverage-cap, --threads, --min-aligned-length (0 = off), --min-aligned-fraction (0 = off),
#   --read-name-regex (accepted for compatibility, ignored)
# --aligner-commandline may contain the literal text $bash_ref_fasta, which is replaced by the
# reference being aligned to. It must write SAM to stdout, e.g. `minimap2 -ax map-ont -Y`.

required=(
  sample_base_name suffix command_mem read_length aligner_commandline
  input_bam mt_cat mt_shifted_cat mt mt_dict mt_shifted_dict mt_interval_list
  coverage_metrics_script
)
optional=(read_name_regex coverage_cap threads min_aligned_length min_aligned_fraction)

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

threads=${threads:-2}
min_aligned_length=${min_aligned_length:-0}
min_aligned_fraction=${min_aligned_fraction:-0}
ALIGNER_VERSION=$(minimap2 --version 2>&1 | head -n1)
export ALIGNER_VERSION

set -o pipefail
set -e

mkdir out

this_sample=out/"${sample_base_name}${suffix}"
this_output_bam_basename=out/"$(basename "${input_bam}" .bam).remap${suffix}"

coverage_cap_args=()
if [[ -n ${coverage_cap:-} ]]; then
  coverage_cap_args=(--COVERAGE_CAP "$coverage_cap")
fi

# Keep only reads whose total aligned length exceeds min_aligned_length and whose aligned fraction
# exceeds min_aligned_fraction (supplementary pieces are summed per read; secondary are ignored).
filter_reads() { # IN.bam OUT.bam
  if [[ $min_aligned_length == 0 && $min_aligned_fraction == 0 ]]; then
    cp "$1" "$2"
    return
  fi
  samtools view -F 0x104 "$1" | awk -v minlen="$min_aligned_length" -v minpc="$min_aligned_fraction" '
    { s = $6; aln = 0; q = 0
      while (match(s, /^[0-9]+[MIDNSHP=X]/)) { l = substr(s, 1, RLENGTH - 1) + 0; op = substr(s, RLENGTH, 1); s = substr(s, RLENGTH + 1)
        if (op == "M" || op == "=" || op == "X") { aln += l; q += l } else if (op == "I" || op == "S") q += l }
      A[$1] += aln; if (q > Q[$1]) Q[$1] = q }
    END { for (r in A) if (A[r] > minlen && Q[r] > 0 && A[r] / Q[r] > minpc) print r }' > keep.txt
  samtools view -b -N keep.txt -o "$2" "$1"
}

# Align reads.fastq to one reference and write <base><tag>.bam (mt contigs only, header reordered to
# the mt dictionary) plus <base><tag>_pre_mt_filt.bam (all contigs), the aligner stderr log and a
# duplicate-metrics placeholder. tag is empty for the unshifted reference and ".shifted" for the shifted one.
align_and_clean() {
  local tag=$1 cat_fasta=$2 dict=$3
  local aligner_cmd=${aligner_commandline//\$bash_ref_fasta/$cat_fasta}

  # shellcheck disable=SC2086  # the command line is intentionally word-split
  $aligner_cmd reads.fastq \
    2> >(tee "${this_output_bam_basename}${tag}.aligner.stderr.log" >&2) |
    samtools view -b -F 0x900 -o "aligned${tag}.bam" -

  filter_reads "aligned${tag}.bam" "filtered${tag}.bam"
  samtools sort -@ "$threads" -o "${this_output_bam_basename}${tag}_pre_mt_filt.bam" "filtered${tag}.bam"
  samtools index "${this_output_bam_basename}${tag}_pre_mt_filt.bam" "${this_output_bam_basename}${tag}_pre_mt_filt.bai"

  # now subset to the mito contigs and use the mt sequence dictionary header
  mapfile -t mt_contigs < <(awk '$1 == "@SQ" {sub("SN:", "", $2); print $2}' "$dict")
  samtools view -h -b -o "mt${tag}.bam" "${this_output_bam_basename}${tag}_pre_mt_filt.bam" "${mt_contigs[@]}"
  samtools view -H "mt${tag}.bam" | grep -v '^@SQ' > "hdr${tag}.sam"
  grep '^@SQ' "$dict" >> "hdr${tag}.sam"
  samtools reheader -P "hdr${tag}.sam" "mt${tag}.bam" > "${this_output_bam_basename}${tag}.bam"
  samtools index "${this_output_bam_basename}${tag}.bam" "${this_output_bam_basename}${tag}.bai"

  printf 'Duplicate marking is not performed for long reads.\n' > "${this_output_bam_basename}${tag}.metrics"
}

# Modified-base tags (MM/ML) and the dorado move table (mv) are kept on the fastq header.
samtools fastq -T MM,ML,mv "${input_bam}" > reads.fastq

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
