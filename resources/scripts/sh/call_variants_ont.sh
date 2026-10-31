#!/bin/bash
# shellcheck disable=SC2154  # the options are assigned dynamically by `declare`
# Long-read (ONT) variant calling on chrM / NUMT regions with Clair3 or ClairS-TO, followed by the
# post-processing that the former GATK scripts (HaplotypeCaller / Mutect2 / FilterMutectCalls) did:
# soft filters (LowAF, LowDP, Homopolymer, blacklisted_site), PASS selection and multi-allelic splitting.
# Shared by MongoHC, MongoNucM2, MongoRunM2InitialFilterSplit and MongoCallMtAndShifted in
# MongoTasks.wdl, so that the output contract of every task stays the same. Outputs go under ./out:
#   <prefix>.raw.vcf (+.idx)          direct caller output (snv + indel merged)
#   <prefix>.raw.vcf.stats            placeholder (Mutect2 stats do not exist for Clair); keeps the layout
#   <prefix>.bamout.bam               empty placeholder (Cromwell dislikes optional outputs)
#   <prefix>.vcf / .filtered.vcf      soft-filtered VCF (same content)
#   <prefix>.pass.vcf                 PASS only
#   <prefix>.pass.split.vcf           PASS only, multi-allelics split
#   <prefix>.splitAndPassOnly.vcf     same as .pass.split.vcf (name kept for haplocheck)
#   <prefix>.passvars.txt             number of PASS records
#   <prefix>.tier<NN>.vcf             optional: .pass.split.vcf re-thresholded at each --tier-thresholds value
#
# Usage: call_variants_ont.sh --out-prefix NAME --caller {clair3|clairs_to} --ref-fasta F --input-bam F ...
# Optional: --interval-list (Picard interval_list or BED), --force-call-vcf, --min-af, --min-dp,
#   --blacklisted-sites (BED), --clair3-model, --clairs-platform, --threads, --extra-args (word-split),
#   --homopolymer-min-len (default 4; 0 disables), --homopolymer-max-af (default 0.5),
#   --tier-thresholds (comma separated AFs, e.g. 0.03,0.05,0.10)

required=(out_prefix caller ref_fasta input_bam)
optional=(
  interval_list force_call_vcf min_af min_dp blacklisted_sites clair3_model clairs_platform
  threads extra_args homopolymer_min_len homopolymer_max_af tier_thresholds
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
case "$caller" in
  clair3 | clairs_to) ;;
  *) echo "ERROR: --caller must be clair3 or clairs_to (got '$caller')" >&2; exit 2 ;;
esac

set -euo pipefail

threads=${threads:-2}
clair3_model=${clair3_model:-r1041_e82_400bps_sup_v520}
clairs_platform=${clairs_platform:-ont_r10_dorado_sup_5khz}
homopolymer_min_len=${homopolymer_min_len:-4}
homopolymer_max_af=${homopolymer_max_af:-0.5}

# Inputs are made absolute and the work happens in a private scratch directory, because the WDL runs
# this script more than once per task (regular + shifted reference).
mkdir -p out
outdir=$(realpath out)
for v in ref_fasta input_bam interval_list force_call_vcf blacklisted_sites; do
  [[ -n ${!v:-} ]] && declare "$v=$(realpath "${!v}")"
done
workdir=$(mktemp -d "$PWD/work.XXXXXX")
trap 'rm -rf -- "$workdir"' EXIT
cd "$workdir"
prefix="${outdir}/${out_prefix}"
touch "${prefix}.bamout.bam"

# --- local, writable copies of the inputs (Cromwell inputs may be read-only symlinks) ---
cp "${ref_fasta}" ref.fa
if [[ -f ${ref_fasta}.fai ]]; then cp "${ref_fasta}.fai" ref.fa.fai; else samtools faidx ref.fa; fi
if [[ -f ${input_bam}.bai ]]; then
  ln -s "${input_bam}" in.bam
  ln -s "${input_bam}.bai" in.bam.bai
elif [[ -f ${input_bam%.bam}.bai ]]; then
  ln -s "${input_bam}" in.bam
  ln -s "${input_bam%.bam}.bai" in.bam.bai
else
  ln -s "${input_bam}" in.bam
  samtools index in.bam
fi

# --- regions: Picard interval_list -> 0-based BED (a BED is used as is) ---
bed_present=false
if [[ -n ${interval_list:-} ]]; then
  if grep -q '^@' "${interval_list}"; then
    awk -F'\t' '!/^@/ {print $1 "\t" $2 - 1 "\t" $3}' "${interval_list}" > regions.bed
  else
    cp "${interval_list}" regions.bed
  fi
  bed_present=true
fi

# --- force calling / genotyping at known sites ---
force_args=()
if [[ -n ${force_call_vcf:-} ]]; then
  echo "Force calling sites from: ${force_call_vcf}"
  bcftools view -Oz -o force.vcf.gz "${force_call_vcf}"
  tabix -f -p vcf force.vcf.gz
  if [[ $caller == clair3 ]]; then force_args=(--vcf_fn=force.vcf.gz); else force_args=(--genotyping_mode_vcf_fn force.vcf.gz); fi
fi

# --- calling ---
# The caller's own AF floor is capped at 0.05: a stricter filter (e.g. 0.95 for nuclear homoplasmic
# calls) is applied afterwards as the LowAF soft filter, not inside the caller.
call_af=""
if [[ -n ${min_af:-} ]]; then
  call_af=$(awk -v a="${min_af}" 'BEGIN {print (a > 0.05 ? 0.05 : a)}')
fi
if [[ $caller == clair3 ]]; then
  # A bare model name is resolved inside the bioconda/docker model directory.
  if [[ $clair3_model == */* ]]; then
    model_path=$clair3_model
  else
    model_path="$(dirname "$(command -v run_clair3.sh)")/models/${clair3_model}"
  fi
  region_args=()
  [[ $bed_present == true ]] && region_args=(--bed_fn=regions.bed)
  af_args=()
  if [[ -n ${call_af:-} ]]; then af_args=(--snp_min_af="${call_af}" --indel_min_af="${call_af}"); fi
  # shellcheck disable=SC2086  # extra_args is intentionally word-split
  run_clair3.sh \
    --bam_fn=in.bam --ref_fn=ref.fa --threads="${threads}" --platform=ont \
    --model_path="${model_path}" --output=clair_out --include_all_ctgs \
    --sample_name="${out_prefix}" --no_phasing_for_fa --remove_intermediate_dir \
    "${region_args[@]}" "${force_args[@]}" "${af_args[@]}" ${extra_args:-}
  cp clair_out/merge_output.vcf.gz raw_in.vcf.gz
  caller_version=$(run_clair3.sh --version 2>&1 | head -n1 || true)
else
  region_args=()
  [[ $bed_present == true ]] && region_args=(--bed_fn regions.bed)
  af_args=()
  if [[ -n ${call_af:-} ]]; then af_args=(--snv_min_af "${call_af}" --indel_min_af "${call_af}"); fi
  # The reference is a single small contig, so one chunk covers it (as MitoHPC2 does).
  chunk=$(awk 'NR == 1 {print $2}' ref.fa.fai)
  # Homoplasmic / near-homoplasmic variants are legitimate here (they feed the self reference),
  # so the somatic-vs-germline tagging modules are switched off.
  # shellcheck disable=SC2086  # extra_args is intentionally word-split
  run_clairs_to \
    --tumor_bam_fn in.bam --ref_fn ref.fa --threads "${threads}" --platform "${clairs_platform}" \
    --output_dir clair_out --include_all_ctgs --sample_name "${out_prefix}" --chunk_size "${chunk}" \
    --disable_nonsomatic_tagging --disable_verdict --remove_intermediate_dir \
    "${region_args[@]}" "${force_args[@]}" "${af_args[@]}" ${extra_args:-}
  for f in clair_out/snv.vcf.gz clair_out/indel.vcf.gz; do
    [[ -f $f ]] || { echo "ERROR: ClairS-TO did not produce $f" >&2; exit 1; }
    tabix -f -p vcf "$f"
  done
  bcftools concat -a -Oz -o raw_in.vcf.gz clair_out/snv.vcf.gz clair_out/indel.vcf.gz
  caller_version=$(run_clairs_to --version 2>&1 | head -n1 || true)
fi
bcftools sort -Oz -o raw_sorted.vcf.gz raw_in.vcf.gz
tabix -f -p vcf raw_sorted.vcf.gz

# Downstream (haplocheck, self reference, merging) needs FORMAT/GT, AF and DP.
header=$(bcftools view -h raw_sorted.vcf.gz)
for tag in GT AF DP; do
  if ! grep -q "^##FORMAT=<ID=${tag}," <<< "$header"; then
    echo "ERROR: caller output lacks FORMAT/${tag}, which downstream steps require" >&2
    exit 1
  fi
done

bcftools view -Ov -o "${prefix}.raw.vcf" raw_sorted.vcf.gz
printf 'caller=%s\nversion=%s\nmin_af=%s\nmin_dp=%s\n' \
  "$caller" "$caller_version" "${min_af:-default}" "${min_dp:-none}" > "${prefix}.raw.vcf.stats"

# --- soft filters (replaces FilterMutectCalls / VariantFiltration) ---
cp "${prefix}.raw.vcf" step.vcf
soft_filter() { # NAME EXPRESSION_TO_EXCLUDE
  bcftools filter -s "$1" -e "$2" -Ov -o step_next.vcf step.vcf
  mv step_next.vcf step.vcf
}
if [[ -n ${min_af:-} ]]; then soft_filter LowAF "MAX(FMT/AF) < ${min_af}"; fi
if [[ -n ${min_dp:-} ]]; then soft_filter LowDP "FMT/DP < ${min_dp}"; fi
if [[ $homopolymer_min_len -gt 0 ]]; then
  # ONT systematically mis-calls the length of homopolymers, so indels inside a run of at least
  # homopolymer_min_len identical bases are dropped unless their AF is high (as in MitoHPC2).
  awk -v n="${homopolymer_min_len}" '
    /^>/ { flush(); name = substr($1, 2); pos = 0; run = 0; prev = ""; next }
    { line = toupper($0); for (i = 1; i <= length(line); i++) { c = substr(line, i, 1); pos++
        if (c == prev) { run++ } else { flush(); prev = c; run = 1; start = pos } } }
    function flush() { if (run >= n && prev != "N") print name "\t" start - 1 "\t" start - 1 + run }
    END { flush() }' ref.fa > hp.bed
  if [[ -s hp.bed ]]; then
    bgzip -c hp.bed > hp.bed.gz
    tabix -f -p bed hp.bed.gz
    printf '##INFO=<ID=HP,Number=0,Type=Flag,Description="Inside a homopolymer run">\n' > hp.hdr
    bgzip -c step.vcf > step.vcf.gz
    bcftools annotate -a hp.bed.gz -c CHROM,FROM,TO -m +HP -h hp.hdr -Ov -o step_next.vcf step.vcf.gz
    mv step_next.vcf step.vcf
    soft_filter Homopolymer "TYPE=\"indel\" && INFO/HP=1 && MAX(FMT/AF) < ${homopolymer_max_af}"
  fi
fi
if [[ -n ${blacklisted_sites:-} ]]; then
  bcftools filter -s blacklisted_site -m + --mask-file "${blacklisted_sites}" -Ov -o step_next.vcf step.vcf
  mv step_next.vcf step.vcf
fi
cp step.vcf "${prefix}.vcf"
cp step.vcf "${prefix}.filtered.vcf"

# --- PASS only, then split multi-allelics (left-aligned, as LeftAlignAndTrimVariants did) ---
bcftools view -f PASS,. -Ov -o "${prefix}.pass.vcf" step.vcf
bcftools view -H "${prefix}.pass.vcf" | wc -l | sed 's/^ *//' > "${prefix}.passvars.txt"
# Splitting happens before the per-allele AF check: on a multi-allelic record the AF floor must apply
# to each ALT allele separately, not to the largest one.
bcftools norm -m -any -f ref.fa -Ov -o split_all.vcf "${prefix}.pass.vcf"
if [[ -n ${min_af:-} ]]; then
  bcftools filter -s LowAF -e "FMT/AF < ${min_af}" -Ov -o split_af.vcf split_all.vcf
else
  cp split_all.vcf split_af.vcf
fi
bcftools view -f PASS,. -Ov -o "${prefix}.pass.split.vcf" split_af.vcf
cp "${prefix}.pass.split.vcf" "${prefix}.splitAndPassOnly.vcf"

# --- optional tiered heteroplasmy thresholds from the same calling pass ---
if [[ -n ${tier_thresholds:-} ]]; then
  IFS=',' read -ra tiers <<< "${tier_thresholds}"
  for t in "${tiers[@]}"; do
    tag=$(awk -v t="$t" 'BEGIN {printf "%02d", t * 100 + 0.5}')
    bcftools view -i "FMT/AF >= ${t}" -Ov -o "${prefix}.tier${tag}.vcf" "${prefix}.pass.split.vcf"
  done
fi

# --- Tribble indices, which the WDL tasks declare as outputs ---
for f in "${prefix}.raw.vcf" "${prefix}.vcf" "${prefix}.filtered.vcf" "${prefix}.pass.vcf" \
  "${prefix}.pass.split.vcf" "${prefix}.splitAndPassOnly.vcf"; do
  gatk IndexFeatureFile -I "$f"
done
