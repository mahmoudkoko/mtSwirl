#!/bin/bash
# shellcheck disable=SC2154  # the required options are assigned dynamically by `declare`
# Build the personalized (self) reference files for one sample.
# Called by the MongoProduceSelfReference task in MongoTasks.wdl.
# All outputs are written under ./out (the WDL collects them from there).
#
# Usage: produce_self_reference.sh --sample-name NAME --input-nuc-vcf F ... (all options required)
# Not using pipefail: several `grep | wc -l` counts below rely on grep returning 1 when nothing matches.
set -eu


required=(
  sample_name suffix n_shift command_mem hail_mem
  input_nuc_vcf input_mt_vcf
  ref_fasta nuc_interval_list
  mt_ref_fasta mt_ref_fasta_index mt_interval_list non_control_region_interval_list
  fa_renaming_script variant_bounds_script check_hom_overlap_script
  hail_force_call_script self_reference_helper_script
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

mt_intervals_basename=$(basename "$mt_interval_list" .interval_list)
noncontrol_basename=$(basename "$non_control_region_interval_list" .interval_list)
selfref_intervals_suffix=.SelfRefLiftover.interval_list

mkdir out

gatk --java-options "-Xmx${command_mem}m" IntervalListTools \
--SORT true \
-I "${mt_interval_list}" \
-O internal_mt.interval_list

gatk --java-options "-Xmx${command_mem}m" IntervalListTools \
--SORT true \
-I "${nuc_interval_list}" \
-O internal_nuc.interval_list

gatk --java-options "-Xmx${command_mem}m" ExtractSequences \
--INTERVAL_LIST internal_mt.interval_list \
-R "${mt_ref_fasta}" \
-O mt_fasta.fasta

Rscript --vanilla "${fa_renaming_script}" mt_fasta.fasta internal_mt.interval_list FALSE mt_fasta_renamed.fasta TRUE
samtools faidx mt_fasta_renamed.fasta

this_nuc_vcf="${input_nuc_vcf}"
this_mt_vcf="${input_mt_vcf}"
this_samp_name="${sample_name}"
this_sample=out/"${sample_name}"

this_mt_chain=out/reference_to_"${this_samp_name}.chain"
this_nuc_chain=out/reference_to_"${this_samp_name}NucOnly.chain"
this_basename="${this_sample}${suffix}"
this_basename_both="${this_sample}andNuc${suffix}"
this_basename_nuc="${this_sample}NucOnly${suffix}"
this_shifted_basename="${this_basename}.shifted_by_${n_shift}_bases"
this_shifted_basename_append="${this_shifted_basename}.cat"

this_filt_vcf="${this_sample}.mtref.homoplasmies.vcf"
this_mt_fasta="${this_basename}.fasta"
this_nuc_mt_fasta="${this_basename_both}.fasta"
this_nuc_only_fasta="${this_basename_nuc}.fasta"
this_fasta_shifted="${this_shifted_basename}.fasta"
this_fasta_cat_shifted="${this_shifted_basename_append}.fasta"
this_chain_shifted="${this_shifted_basename}.shift_back_${n_shift}_bases.chain"
this_chain_fwd_shifted="${this_shifted_basename}.shift_fwd_${n_shift}_bases.chain"
this_mt_intervals=out/"${mt_intervals_basename}.${this_samp_name}${selfref_intervals_suffix}"
this_nonctrl_interval=out/"${noncontrol_basename}.${this_samp_name}${selfref_intervals_suffix}"
this_ctrl_interval=out/control_region_shifted.chrM."${this_samp_name}${selfref_intervals_suffix}"

this_vcf_bgz="${this_basename}.reversed.selfRef.homoplasmies.vcf.bgz"
this_vcf_filters_bgz="${this_basename}.reversed.withfilters.selfRef.homoplasmies.vcf.bgz"
this_vcf_shifted_bgz="${this_basename}.reversed.selfRef.shifted.homoplasmies.vcf.bgz"

bgzip -c "${this_mt_vcf}" > input_mt_vcf.gz && tabix input_mt_vcf.gz
bcftools view -Oz -i 'FORMAT/AF>0.95' input_mt_vcf.gz > "${this_filt_vcf}.gz"
tabix "${this_filt_vcf}.gz"
bgzip -cd "${this_filt_vcf}.gz" > "${this_filt_vcf}"

Rscript --vanilla "${check_hom_overlap_script}" "${this_sample}" "${this_filt_vcf}" "${this_nuc_vcf}" FALSE
bgzip -c "overlapping_variants_to_remove.vcf" > mt_overlaps_rm.vcf.gz && tabix mt_overlaps_rm.vcf.gz
bcftools isec "${this_filt_vcf}.gz" mt_overlaps_rm.vcf.gz -p output_isec
private_to_rem=$(cat ./output_isec/0001.vcf | grep ^chr | wc -l | sed 's/^ *//g')
export private_to_rem
if [ "$private_to_rem" -ne 0 ]; then
  echo "ERROR: There should not be any variants private to the VCF for removal."
  exit 1;
fi
shared_by_both=$(cat ./output_isec/0002.vcf | grep ^chr | wc -l | sed 's/^ *//g')
export shared_by_both
for_removal=$(cat overlapping_variants_to_remove.vcf | grep ^chr | wc -l | sed 's/^ *//g')
export for_removal
if [ "$shared_by_both" -ne "$for_removal" ]; then
  echo "ERROR: The number of variants shared by both rejected VCF and MT VCF should be the number of variants for removal."
  exit 1;
fi
only_undupe=$(cat ./output_isec/0000.vcf | grep ^chr | wc -l | sed 's/^ *//g')
export only_undupe
original_nrow=$(cat "${this_filt_vcf}" | grep ^chr | wc -l | sed 's/^ *//g')
export original_nrow
if [ "$((original_nrow-for_removal))" -ne "$only_undupe" ]; then
  echo "ERROR: New VCF should be smaller than old VCF by the exact number of records for removal."
  exit 1;
fi
Rscript --vanilla "${check_hom_overlap_script}" "temp" "output_isec/0000.vcf" "${this_nuc_vcf}" TRUE
rm "${this_filt_vcf}.gz" "${this_filt_vcf}" "${this_filt_vcf}.gz.tbi" mt_overlaps_rm.vcf.gz mt_overlaps_rm.vcf.gz.tbi overlapping_variants_to_remove.vcf
cp "output_isec/0000.vcf" "${this_filt_vcf}"
bgzip -c "${this_filt_vcf}" > "${this_filt_vcf}".gz && tabix "${this_filt_vcf}".gz
rm -rf output_isec

# Select the high-AF sample's ALT allele rather than relying on bcftools defaults.
bcftools consensus -s "${this_samp_name}" -H A -f mt_fasta_renamed.fasta -o mt_fasta_lifted.fasta -c "${this_mt_chain}" "${this_filt_vcf}.gz"
Rscript --vanilla "${fa_renaming_script}" mt_fasta_lifted.fasta internal_mt.interval_list TRUE "${this_mt_fasta}" TRUE

gatk --java-options "-Xmx${command_mem}m" CreateSequenceDictionary \
-R "${this_mt_fasta}" \
-O "${this_basename}.dict"

samtools faidx "${this_mt_fasta}"

Rscript --vanilla "${variant_bounds_script}" "${this_nuc_vcf}" internal_nuc.interval_list rejected_nuc.vcf
bgzip -c "${this_nuc_vcf}" > input_nuc_vcf.gz && tabix input_nuc_vcf.gz
bgzip -c rejected_nuc.vcf > rejected_nuc.vcf.gz && tabix rejected_nuc.vcf.gz
bcftools isec -p intersected_vcfs -Ov input_nuc_vcf.gz rejected_nuc.vcf.gz
private_to_rej=$(cat ./intersected_vcfs/0001.vcf | grep ^chr | wc -l | sed 's/^ *//g')
export private_to_rej
if [ "$private_to_rej" -ne 0 ]; then
  echo "ERROR: There should not be any variants private to the rejected VCF."
  exit 1;
fi
cat ./intersected_vcfs/0002.vcf | grep ^chr | wc -l | sed 's/^ *//g' > "${this_sample}.nuc.removed.txt"

gatk --java-options "-Xmx${command_mem}m" ExtractSequences \
--INTERVAL_LIST internal_nuc.interval_list \
-R "${ref_fasta}" \
-O nuc_fasta.fasta

Rscript --vanilla "${fa_renaming_script}" nuc_fasta.fasta internal_nuc.interval_list FALSE nuc_fasta_renamed.fasta FALSE
samtools faidx nuc_fasta_renamed.fasta
bgzip -c ./intersected_vcfs/0000.vcf > filtered_nuc_vcf.gz && tabix filtered_nuc_vcf.gz
bcftools consensus -s "${this_samp_name}" -H A -f nuc_fasta_renamed.fasta -o nuc_fasta_lifted.fasta -c "${this_nuc_chain}" filtered_nuc_vcf.gz
Rscript --vanilla "${fa_renaming_script}" nuc_fasta_lifted.fasta internal_nuc.interval_list TRUE "${this_nuc_only_fasta}" FALSE
cat "${this_nuc_only_fasta}" "${this_mt_fasta}" > "${this_nuc_mt_fasta}"

gatk --java-options "-Xmx${command_mem}m" CreateSequenceDictionary \
-R "${this_nuc_mt_fasta}" \
-O "${this_basename_both}.dict"

samtools faidx "${this_nuc_mt_fasta}"

gatk --java-options "-Xmx${command_mem}m" CreateSequenceDictionary \
-R "${this_nuc_only_fasta}" \
-O "${this_basename_nuc}.dict"

samtools faidx "${this_nuc_only_fasta}"

echo "Now shifting the reference..."
Rscript --vanilla "${self_reference_helper_script}" shift-reference \
  "${this_mt_fasta}" "${n_shift}" \
  "${this_fasta_shifted}" \
  "${this_chain_shifted}" \
  "${this_chain_fwd_shifted}"

cat "${this_fasta_shifted}" "${this_nuc_only_fasta}" > "${this_fasta_cat_shifted}"

gatk --java-options "-Xmx${command_mem}m" CreateSequenceDictionary \
-R "${this_fasta_cat_shifted}" \
-O "${this_shifted_basename_append}.dict"

samtools faidx "${this_fasta_cat_shifted}"

gatk --java-options "-Xmx${command_mem}m" CreateSequenceDictionary \
-R "${this_fasta_shifted}" \
-O "${this_shifted_basename}.dict"

samtools faidx "${this_fasta_shifted}"

echo "Now adjusting the MT interval list..."
Rscript --vanilla "${self_reference_helper_script}" adjust-mt-intervals \
  "${mt_interval_list}" \
  "${this_basename}.dict" \
  "${this_mt_intervals}"

echo "Now shifting the noncontrol region..."

gatk --java-options "-Xmx${command_mem}m" LiftOverIntervalList \
-I "${non_control_region_interval_list}" \
-O "${this_nonctrl_interval}" \
-SD "${this_basename}.dict" \
--CHAIN "${this_mt_chain}"

echo "Now shifting the control region..."
Rscript --vanilla "${self_reference_helper_script}" shift-control-region \
  "${this_nonctrl_interval}" \
  "${this_shifted_basename}.dict" \
  "${n_shift}" \
  "${this_ctrl_interval}"

echo "Now making force call variants..."
export PYSPARK_SUBMIT_ARGS="--driver-memory ${hail_mem}g pyspark-shell"
python "${hail_force_call_script}" \
  --self-fasta "${this_mt_fasta}" \
  --self-fai "${this_mt_fasta}.fai" \
  --reference-fasta "${mt_ref_fasta}" \
  --reference-fai "${mt_ref_fasta_index}" \
  --reference-to-self-chain "${this_mt_chain}" \
  --shifted-self-fasta "${this_fasta_shifted}" \
  --shifted-self-fai "${this_fasta_shifted}.fai" \
  --forward-shift-chain "${this_chain_fwd_shifted}" \
  --filtered-vcf "${this_filt_vcf}" \
  --output-vcf "${this_vcf_bgz}" \
  --output-filter-vcf "${this_vcf_filters_bgz}" \
  --output-shifted-vcf "${this_vcf_shifted_bgz}"
