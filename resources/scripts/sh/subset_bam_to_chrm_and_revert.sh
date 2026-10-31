#!/bin/bash
# shellcheck disable=SC2154  # the options are assigned dynamically by `declare`
# Long-read (ONT) version. Subset a WGS BAM/CRAM to the mt (and optionally NUMT) intervals, keep
# primary alignments only, drop reads that are mostly NUMT / chimeric (aligned length and aligned
# fraction filter, as in MitoHPC2 filter.lr.sh), optionally cap the number of reads, revert to an
# unmapped BAM, and collect WGS metrics and mean/median coverage. There is no mate rescue and no
# duplicate marking for long reads; the sorted .proc.bam is the filtered subset itself. Called by MongoSubsetBamToChrMAndRevert
# in MongoTasks.wdl. All outputs are written under ./out (the WDL collects them from there).
#
# Usage: subset_bam_to_chrm_and_revert.sh --sample-name NAME ... (see `required`)
# Optional: --ref-fasta, --mt-interval-list, --nuc-interval-list, --contig-name,
#   --requester-pays-project, --coverage-cap, --min-aligned-length (default 0 = off),
#   --min-aligned-fraction (default 0 = off), --max-reads (default 0 = no cap),
#   --read-name-regex and --printreads-extra-args (accepted for compatibility, ignored)

required=(
  sample_name command_mem read_length input_bam input_bai
  skip_restore_hardclips coverage_metrics_script
)
optional=(
  ref_fasta mt_interval_list nuc_interval_list contig_name requester_pays_project
  read_name_regex coverage_cap printreads_extra_args
  min_aligned_length min_aligned_fraction max_reads
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

set -euo pipefail

mkdir out

this_sample=out/"${sample_name}"
java_opts="-Xmx${command_mem}m"
min_aligned_length=${min_aligned_length:-0}
min_aligned_fraction=${min_aligned_fraction:-0}
max_reads=${max_reads:-0}

# --- regions: mt + NUMT interval lists (Picard) merged into one 0-based BED, plus an optional contig ---
: > regions.raw.bed
for il in "${mt_interval_list:-}" "${nuc_interval_list:-}"; do
  [[ -n $il ]] && awk -F'\t' '!/^@/ {print $1 "\t" $2 - 1 "\t" $3}' "$il" >> regions.raw.bed
done
[[ -n ${contig_name:-} ]] && printf '%s\t0\t%s\n' "$contig_name" 2147483647 >> regions.raw.bed
region_args=()
if [[ -s regions.raw.bed ]]; then
  sort -k1,1 -k2,2n regions.raw.bed | bedtools merge -i - > regions.bed
  region_args=(-L regions.bed)
fi
ref_args=()
[[ -n ${ref_fasta:-} ]] && ref_args=(-T "$ref_fasta")
[[ -n ${requester_pays_project:-} ]] && export GCS_REQUESTER_PAYS_PROJECT="$requester_pays_project"

# Primary alignments only (drops secondary 0x100 and supplementary 0x800). Local inputs are linked
# next to each other so that the index is found; remote inputs are passed with an explicit index.
if [[ -f $input_bam && -f $input_bai ]]; then
  case "$input_bam" in
    *.cram) ln -s "$(realpath "$input_bam")" in.cram; ln -s "$(realpath "$input_bai")" in.cram.crai; in=in.cram ;;
    *) ln -s "$(realpath "$input_bam")" in.bam; ln -s "$(realpath "$input_bai")" in.bam.bai; in=in.bam ;;
  esac
else
  in="${input_bam}##idx##${input_bai}"
fi
samtools view -b -F 0x900 "${ref_args[@]}" "${region_args[@]}" -o subset.bam "$in"

# Aligned-length / aligned-fraction read filter. A read is kept if its aligned bases exceed
# min_aligned_length (set above the longest NUMT) and make up more than min_aligned_fraction of
# the read, which removes NUMT-derived and chimeric reads.
if [[ $min_aligned_length != 0 || $min_aligned_fraction != 0 ]]; then
  samtools view subset.bam | awk -v minlen="$min_aligned_length" -v minpc="$min_aligned_fraction" '
    { s = $6; aln = 0; q = 0
      while (match(s, /^[0-9]+[MIDNSHP=X]/)) { l = substr(s, 1, RLENGTH - 1) + 0; op = substr(s, RLENGTH, 1); s = substr(s, RLENGTH + 1)
        if (op == "M" || op == "=" || op == "X") { aln += l; q += l } else if (op == "I" || op == "S") q += l }
      A[$1] += aln; if (q > Q[$1]) Q[$1] = q }
    END { for (r in A) if (A[r] > minlen && Q[r] > 0 && A[r] / Q[r] > minpc) print r }' > keep.txt
  samtools view -b -N keep.txt -o filtered.bam subset.bam
else
  mv subset.bam filtered.bam
fi

# Read cap: random downsampling to roughly max_reads reads (deterministic seed).
n_reads=$(samtools view -c filtered.bam)
if [[ $max_reads -gt 0 && $n_reads -gt $max_reads ]]; then
  frac=$(awk -v m="$max_reads" -v n="$n_reads" 'BEGIN {printf "%.6f", m / n}')
  echo "Downsampling ${n_reads} reads to ~${max_reads} (fraction ${frac})"
  samtools view -b -s "42${frac#0}" -o capped.bam filtered.bam
  mv capped.bam filtered.bam
fi

samtools sort -@ 2 -o "${this_sample}.bam" filtered.bam
samtools index "${this_sample}.bam"
echo 0 > "${this_sample}.ct_failed.txt"

echo "Now removing mapping..."
samtools reset -o "${this_sample}.unmap.bam" "${this_sample}.bam"

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

# No duplicate marking for long reads: keep the output contract with a placeholder metrics file.
printf 'Duplicate marking is not performed for long reads.\n' > "${this_sample}.duplicate.metrics"
cp "${this_sample}.bam" "${this_sample}.proc.bam"
samtools index "${this_sample}.proc.bam" "${this_sample}.proc.bai"
