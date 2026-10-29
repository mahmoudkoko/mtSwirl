#!/bin/bash
# shellcheck disable=SC2154  # the options are assigned dynamically by `declare`
# Estimate mt contamination with Haplochecker and write one small text file per value
# (contamination flag, major/minor haplogroup, mean het levels). Called by GetContamination in
# MongoTasks.wdl. All outputs are written under ./out (the WDL collects them from there).
#
# Usage: get_contamination.sh --sample-name NAME --command-mem MB --mean-coverage N \
#          --input-vcf VCF --haplocheck-cli JAR

required=(sample_name command_mem mean_coverage input_vcf haplocheck_cli)

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

# Not using pipefail: the variant count below relies on grep returning 1 when nothing matches.
set -e

mkdir out
this_basename=out/"${sample_name}"

this_vcf_nvar=$(grep -c '^chrM' "${input_vcf}" || true)
echo "${sample_name} has VCF with ${this_vcf_nvar} variants for contamination."

# Haplochecker reads every VCF in the input directory and writes ./output
java -Xmx"${command_mem}"m -jar "${haplocheck_cli}" "$(dirname "${input_vcf}")"

sed 's/\"//g' output > output-noquotes
cp 'output-noquotes' "${this_basename}_output_noquotes"

grep "SampleID" output-noquotes > headers
FORMAT_ERROR="Bad contamination file format"
check_header() {
  if [ "$(awk -v c="$1" '{print $c}' headers)" != "$2" ]; then
    echo "$FORMAT_ERROR"
    exit 1
  fi
}
check_header 2 Contamination
check_header 6 HgMajor
check_header 8 HgMinor
check_header 14 MeanHetLevelMajor
check_header 15 MeanHetLevelMinor

if [ "${mean_coverage}" -gt 0 ] && [ "${this_vcf_nvar}" -gt 0 ]; then
  grep -v "SampleID" output-noquotes > output-data
  awk -F "\t" '{print $2}' output-data > "${this_basename}.contamination.txt"
  awk -F "\t" '{print $6}' output-data > "${this_basename}.major_hg.txt"
  awk -F "\t" '{print $8}' output-data > "${this_basename}.minor_hg.txt"
  awk -F "\t" '{print $14}' output-data > "${this_basename}.mean_het_major.txt"
  awk -F "\t" '{print $15}' output-data > "${this_basename}.mean_het_minor.txt"
else
  if [ "${mean_coverage}" -eq 0 ]; then
    echo "Sample ${sample_name} has a mean coverage of 0."
  elif [ "${this_vcf_nvar}" -eq 0 ]; then
    echo "Sample ${sample_name} has no variants."
  else
    echo "Unsupported contamination exception for ${sample_name}"
    exit 1
  fi

  echo "NO" > "${this_basename}.contamination.txt"
  echo "NONE" > "${this_basename}.major_hg.txt"
  echo "NONE" > "${this_basename}.minor_hg.txt"
  echo "0.000" > "${this_basename}.mean_het_major.txt"
  echo "0.000" > "${this_basename}.mean_het_minor.txt"
fi
