#!/bin/bash
# shellcheck disable=SC2154  # the options are assigned dynamically by `declare`
# Lift the second-round self-reference variants back to GRCh38 (LiftoverVcf, then the Hail
# rescue of rejected records), merge them, and build the per-base coverage table.
# Called by MongoLiftoverVCFAndGetCoverage in MongoTasks.wdl.
# All outputs are written under ./out (the WDL collects them from there).
#
# Usage: liftover_vcf_and_coverage.sh --sample-name NAME --self-suffix S ... (see `required`)
# Not using pipefail: several `grep | wc -l` counts below rely on grep returning 1 when nothing matches.
set -e

required=(
  sample_name self_suffix hail_mem
  original_filtered_vcf new_self_ref_vcf reversed_hom_ref_vcf
  mt_self mt_self_index mt_self_shifted
  chain_self_to_ref chain_ref_to_self
  input_bam_regular_ref input_bam_shifted_ref
  self_control_region_shifted_reference_interval_list self_non_control_region_interval_list
  hail_liftover coverage_helpers_script
  ref_fasta ref_fasta_index
)

while [[ $# -gt 0 ]]; do
  key=${1#--}
  key=${key//-/_}
  if [[ $1 != --* || $# -lt 2 || ! " ${required[*]} " == *" $key "* ]]; then
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

mkdir out

this_sample_name="${sample_name}"
this_self_ref_vcf="${new_self_ref_vcf}"
this_ref_filtered_vcf="${original_filtered_vcf}"
this_rev_hom_ref_vcf="${reversed_hom_ref_vcf}"
this_self_to_ref_chain="${chain_self_to_ref}"
this_ref_to_self_chain="${chain_ref_to_self}"
this_self_fasta="${mt_self}"
this_self_fai="${mt_self_index}"
this_self_shifted_fasta="${mt_self_shifted}"
this_self_control="${self_control_region_shifted_reference_interval_list}"
this_self_non_control="${self_non_control_region_interval_list}"
this_self_bam="${input_bam_regular_ref}"
this_self_shifted_bam="${input_bam_shifted_ref}"

this_sample=out/"${this_sample_name}"
this_basename="${this_sample}${self_suffix}.split"
this_logging="${this_basename}_fix_liftover.log"

bgzip -c "${this_self_ref_vcf}" > "${this_self_ref_vcf}.bgz" && tabix "${this_self_ref_vcf}.bgz"
tabix "${this_rev_hom_ref_vcf}"
bcftools isec -p intersected_vcfs -Ov "${this_self_ref_vcf}.bgz" "${this_rev_hom_ref_vcf}"

# there should be no records private to reversed hom ref VCF
private_to_rev_hom_ref=$(cat ./intersected_vcfs/0001.vcf | grep ^chrM | wc -l | sed 's/^ *//g')
export private_to_rev_hom_ref
if [ "$private_to_rev_hom_ref" -ne 0 ]; then
  echo "ERROR: There should not be any variants private to the reversed hom ref VCF."
  exit 1;
fi

gatk LiftoverVcf \
  -I ./intersected_vcfs/0000.vcf \
  -O "${this_basename}.selfToRef.pre.vcf" \
  -R "${ref_fasta}" \
  --CHAIN "${this_self_to_ref_chain}" \
  --REJECT "${this_basename}.selfToRef.rejected.pre.vcf"

gatk MergeVcfs \
  -I "${this_basename}.selfToRef.rejected.pre.vcf" \
  -I ./intersected_vcfs/0002.vcf \
  -O "${this_basename}.selfToRef.rejected.vcf"

sed -e 's/^chr//' "${this_ref_filtered_vcf}" \
  | awk '{OFS="\t"; if (!/^#/){print $1,$2-(length($4)>1 ? 0 : 1),$2-1+length($4),$4"/"$5,"+",length($4),length($5)}}' > filtered.bed
sed -e 's/^chr//' "${this_basename}.selfToRef.pre.vcf" \
  | awk '{OFS="\t"; if (!/^#/){print $1,$2-(length($4)>1 ? 0 : 1),$2-1+length($4),$4"/"$5,"+",length($4),length($5)}}' > success.bed
n_ref_pass_thru=$(bedtools intersect -a filtered.bed -b success.bed | awk '{OFS="\t"; if (($7 == 1)) {print}}' | wc -l | sed 's/^ *//g')
export n_ref_pass_thru
echo "$n_ref_pass_thru" > "${this_sample}_n_ref_pass_thru.txt"

if [ "$n_ref_pass_thru" -ne 0 ]; then
  echo "ERROR: All variants changed in the self-reference, and all sites within self-reference insertions (excluding the first base), should have failed LiftoverVcf, which is not the case here."
  exit 1;
fi

n_filtered=$(cat "${this_ref_filtered_vcf}" | grep ^chrM | wc -l | sed 's/^ *//g')
export n_filtered
n_original=$(cat "${this_self_ref_vcf}" | grep ^chrM | wc -l | sed 's/^ *//g')
export n_original
n_pass=$(cat "${this_basename}.selfToRef.pre.vcf" | grep ^chrM | wc -l | sed 's/^ *//g')
export n_pass
n_rejected=$(cat "${this_basename}.selfToRef.rejected.vcf" | grep ^chrM | wc -l | sed 's/^ *//g')
export n_rejected
n_pass_rejected=$((n_pass + n_rejected))
export n_pass_rejected
echo "$n_pass" > "${this_sample}_n_pass.txt"

if [ "$n_filtered" -gt "$n_rejected" ]; then
  echo "ERROR: the number of sites changed in the self-reference should be less than or equal to the number of sites rejected on initial LiftoverVcf run."
  exit 1;
fi

if [ "$n_pass_rejected" -ne "$n_original" ]; then
  echo "ERROR: Records appear to have been lost. The sum of records in passing and rejected VCFs should be the same as the number in the original VCF."
  exit 1;
fi

echo "$n_filtered sites are changed in self-reference. $n_original sites were found variant after second-round Mutect. Of these, $n_pass passed first-round Liftover and $n_rejected failed and are being piped to Hail pipeline for rescue."

# now run hail script to fix the rejects
export PYSPARK_SUBMIT_ARGS="--driver-memory ${hail_mem}g pyspark-shell"
python "${hail_liftover}" \
--vcf-file "${this_basename}.selfToRef.rejected.vcf" \
--success-vcf-file "${this_basename}.selfToRef.pre.vcf" \
--self-homoplasmies "${this_rev_hom_ref_vcf}" \
--individual-name "${this_sample_name}" \
--self-to-ref-chain "${this_self_to_ref_chain}" \
--ref-to-self-chain "${this_ref_to_self_chain}" \
--self-fasta "${this_self_fasta}" \
--self-fai "${this_self_fai}" \
--reference-fasta "${ref_fasta}" \
--reference-fai "${ref_fasta_index}" \
--output-prefix "${this_basename}.round2liftover" \
--export-homoplasmic-deletions-coverage \
--output-txt-for-wdl \
--logging "${this_logging}"

bgzip -cd "${this_basename}.round2liftover.rejected.vcf.bgz" > "${this_basename}.round2liftover.rejected.vcf"
bgzip -cd "${this_basename}.round2liftover.fixed.vcf.bgz" > "${this_basename}.round2liftover.fixed.vcf"
bgzip -cd "${this_basename}.round2liftover.updated_success.vcf.bgz" > "${this_basename}.round2liftover.updated_success.vcf"

gatk MergeVcfs \
  -I "${this_basename}.round2liftover.updated_success.vcf" \
  -I "${this_basename}.round2liftover.fixed.vcf" \
  -O "${this_basename}.selfToRef.final.vcf"

n_final_pass=$(cat "${this_basename}.selfToRef.final.vcf" | grep ^chrM | wc -l | sed 's/^ *//g')
export n_final_pass
echo "$n_final_pass" > "${this_sample}_n_final_pass.txt"

echo "Now producing coverage file..."
gatk CollectHsMetrics \
  -I "${this_self_bam}" \
  -R "${this_self_fasta}" \
  --PER_BASE_COVERAGE non_control_region.tsv \
  -O non_control_region.metrics \
  -TI "${this_self_non_control}" \
  -BI "${this_self_non_control}" \
  --COVERAGE_CAP 20000 \
  --SAMPLE_SIZE 1 \
  --VALIDATION_STRINGENCY LENIENT

gatk CollectHsMetrics \
  -I "${this_self_shifted_bam}" \
  -R "${this_self_shifted_fasta}" \
  --PER_BASE_COVERAGE control_region_shifted.tsv \
  -O control_region_shifted.metrics \
  -TI "${this_self_control}" \
  -BI "${this_self_control}" \
  --COVERAGE_CAP 20000 \
  --SAMPLE_SIZE 1 \
  --VALIDATION_STRINGENCY LENIENT

Rscript --vanilla "${coverage_helpers_script}" stitch-self-coverage \
  "${this_self_fasta}" \
  control_region_shifted.tsv \
  non_control_region.tsv \
  "${this_basename}.per_base_coverage.tsv"

echo "Now outputting a final table with integers..."
paste -d "\t" "${this_basename}.round2liftover.all_int_outputs.txt" <(printf "n_liftover_changed_selfref_and_passed\n$(cat "${this_sample}_n_ref_pass_thru.txt")\n") <(printf "n_liftover_r1_pass\n$(cat "${this_sample}_n_pass.txt")\n") <(printf "n_liftover_r2_pass\n$(cat "${this_sample}_n_final_pass.txt")\n") > "${this_basename}.round2liftover.all_int_outputs.final.txt"
