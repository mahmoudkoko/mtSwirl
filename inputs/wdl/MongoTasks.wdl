version 1.0
# This is a stupid file containing giant tasks for use with pipeline platforms
# that cannot handle a lot of parallel complexity. This is obviously hard to maintain.
# These tasks are serial combinations of more modular tasks that were located in each individual WDL.


#################################
# MongoSubsetBamToChrMAndRevert #
#################################

task MongoSubsetBamToChrMAndRevert {
  input {
    File input_bam
    File? input_bai
    String sample_name

    File? mt_interval_list
    File? nuc_interval_list
    String? contig_name
    String? requester_pays_project
    File? ref_fasta
    File? ref_fasta_index
    File? ref_dict
    File coverage_metrics_script
    File subset_script

    Boolean skip_restore_hardclips
    String? read_name_regex
    Int? read_length
    Int? coverage_cap

    String? printreads_extra_args
    Int min_aligned_length = 0
    Float min_aligned_fraction = 0
    Int max_reads = 0
    Int mem
  }
  String appended_bam = input_bam + '.crai'
  Int read_length_for_optimization = select_first([read_length, 151])
  Int command_mem = (mem * 1000) - 500

  meta {
    description: "Subsets a whole genome long-read bam to mt/NUMT primary alignments, filters NUMT-like reads (aligned length / fraction), caps the read count and reverts to an unmapped bam"
  }
  parameter_meta {
    ref_fasta: "Reference is only required for cram input. If it is provided ref_fasta_index and ref_dict are also required."
    input_bam: {
      localization_optional: true
    }
    input_bai: {
      localization_optional: true
    }
  }
  command <<<
    bash ~{subset_script} \
      --sample-name "~{sample_name}" \
      --command-mem ~{command_mem} \
      --read-length ~{read_length_for_optimization} \
      --input-bam "~{input_bam}" \
      --input-bai "~{select_first([input_bai, appended_bam])}" \
      --skip-restore-hardclips ~{skip_restore_hardclips} \
      --coverage-metrics-script "~{coverage_metrics_script}" \
      ~{"--ref-fasta " + ref_fasta} \
      ~{"--mt-interval-list " + mt_interval_list} \
      ~{"--nuc-interval-list " + nuc_interval_list} \
      ~{"--contig-name " + contig_name} \
      ~{"--requester-pays-project " + requester_pays_project} \
      ~{"--read-name-regex " + read_name_regex} \
      ~{"--coverage-cap " + coverage_cap} \
      --min-aligned-length ~{min_aligned_length} \
      --min-aligned-fraction ~{min_aligned_fraction} \
      --max-reads ~{max_reads}
  >>>
  output {
    File output_bam = "out/~{sample_name}.proc.bam"
    File output_bai = "out/~{sample_name}.proc.bai"
    File unmapped_bam = "out/~{sample_name}.unmap.bam"
    File duplicate_metrics = "out/~{sample_name}.duplicate.metrics"
    Int reads_dropped = read_int("out/~{sample_name}.ct_failed.txt")
    Int mean_coverage = read_int("out/~{sample_name}.mean_coverage.txt")
  }
}


###########
# MongoHC #
###########

task MongoHC {
  input {
    File ref_fasta
    File ref_fai
    File ref_dict

    File input_bam
    File input_bai

    String sample_name
    String suffix = ""

    Int max_reads_per_alignment_start = 75
    String? hc_extra_args
    Boolean make_bamout = false

    File? nuc_interval_list
    File? force_call_vcf
    File? force_call_vcf_index
    Float? contamination

    Int haplotype_caller_nucdna_dp_lower_bound
    File hc_script

    String mt_caller = "clairs_to"
    String clair3_model = "r1041_e82_400bps_sup_v520"
    String clairs_platform = "ont_r10_dorado_sup_5khz"
    Int n_cpu = 2
    Int mem
  }

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int command_mem = (mem * 1000) - 500

  command <<<
    bash ~{hc_script} \
      --out-prefix "out/~{sample_name}~{suffix}" \
      --ref-fasta "~{ref_fasta}" \
      --input-bam "~{input_bam}" \
      --min-dp ~{haplotype_caller_nucdna_dp_lower_bound} \
      --caller "~{mt_caller}" \
      --clair3-model "~{clair3_model}" \
      --clairs-platform "~{clairs_platform}" \
      --threads ~{n_cpu} \
      ~{"--interval-list " + nuc_interval_list} \
      ~{"--force-call-vcf " + force_call_vcf} \
      ~{"--extra-args '" + hc_extra_args + "'"}
  >>>

  output {
    File raw_vcf = "out/~{sample_name}~{suffix}.raw.vcf"
    File raw_vcf_idx = "out/~{sample_name}~{suffix}.raw.vcf.idx"
    File output_bamOut = "out/~{sample_name}~{suffix}.bamout.bam"
    File filtered_vcf = "out/~{sample_name}~{suffix}.vcf"
    File filtered_vcf_idx = "out/~{sample_name}~{suffix}.vcf.idx"
    File full_pass_vcf = "out/~{sample_name}~{suffix}.pass.vcf"
    File full_pass_vcf_index = "out/~{sample_name}~{suffix}.pass.vcf.idx"
    Int post_filt_vars = read_int('out/~{sample_name}~{suffix}.passvars.txt')
    File split_vcf = "out/~{sample_name}~{suffix}.pass.split.vcf"
    File split_vcf_index = "out/~{sample_name}~{suffix}.pass.split.vcf.idx"
  }
}


##############
# MongoNucM2 #
##############

task MongoNucM2 {
  input {
    File ref_fasta
    File ref_fai
    File ref_dict
    File input_bam
    File input_bai

    String sample_name
    String suffix = ""

    Int max_reads_per_alignment_start = 75
    String? m2_extra_args
    Boolean make_bamout = false
    File? mt_interval_list
    String? m2_extra_filtering_args
    Int max_alt_allele_count
    Float? vaf_filter_threshold
    Float? f_score_beta
    Float nuc_contamination_estimate
    File m2_script
    File? blacklisted_sites
    File? blacklisted_sites_index
    String mt_caller = "clairs_to"
    String clair3_model = "r1041_e82_400bps_sup_v520"
    String clairs_platform = "ont_r10_dorado_sup_5khz"
    Int n_cpu = 2
    Int mem
  }

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int command_mem = (mem * 1000) - 500


  meta {
    description: "Clair3 / ClairS-TO for calling Snps and Indels (nuclear homoplasmic calls)"
  }
  parameter_meta {
    input_bam: "Aligned Bam"
  }
  command <<<
    bash ~{m2_script} \
      --out-prefix "out/~{sample_name}~{suffix}" \
      --ref-fasta "~{ref_fasta}" \
      --input-bam "~{input_bam}" \
      --caller "~{mt_caller}" \
      --clair3-model "~{clair3_model}" \
      --clairs-platform "~{clairs_platform}" \
      --threads ~{n_cpu} \
      ~{"--interval-list " + mt_interval_list} \
      ~{"--min-af " + vaf_filter_threshold} \
      ~{"--blacklisted-sites " + blacklisted_sites} \
      ~{"--extra-args '" + m2_extra_args + "'"}
  >>>
  output {
    File raw_vcf = "out/~{sample_name}~{suffix}.raw.vcf"
    File raw_vcf_idx = "out/~{sample_name}~{suffix}.raw.vcf.idx"
    File stats = "out/~{sample_name}~{suffix}.raw.vcf.stats"
    File output_bamOut = "out/~{sample_name}~{suffix}.bamout.bam"

    File filtered_vcf = "out/~{sample_name}~{suffix}.vcf"
    File filtered_vcf_idx = "out/~{sample_name}~{suffix}.vcf.idx"

    File full_pass_vcf = "out/~{sample_name}~{suffix}.pass.vcf"
    File full_pass_vcf_index = "out/~{sample_name}~{suffix}.pass.vcf.idx"
    Int post_filt_vars = read_int('out/~{sample_name}~{suffix}.passvars.txt')

    File split_vcf = "out/~{sample_name}~{suffix}.pass.split.vcf"
    File split_vcf_index = "out/~{sample_name}~{suffix}.pass.split.vcf.idx"
  }
}


################################
# MongoRunM2InitialFilterSplit #
################################

task MongoRunM2InitialFilterSplit {
  # intended for use prior to contamination (e.g., assumes mtDNA, no contamination)
  input {
    Int mean_coverage
    String sample_name
    File input_bam
    File input_bai
    Float nuc_contamination_estimate
    File m2_script
    String suffix

    File ref_fasta
    File ref_fai
    File ref_dict
    Int max_reads_per_alignment_start = 75
    String? m2_extra_args
    Boolean make_bamout = false
    File? mt_interval_list
    String? m2_extra_filtering_args
    Int max_alt_allele_count
    Float? vaf_filter_threshold
    Float? f_score_beta

    File? blacklisted_sites
    File? blacklisted_sites_index
    String mt_caller = "clairs_to"
    String clair3_model = "r1041_e82_400bps_sup_v520"
    String clairs_platform = "ont_r10_dorado_sup_5khz"
    Int n_cpu = 2
    Int mem
  }

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int command_mem = (mem * 1000) - 500


  meta {
    description: "Clair3 / ClairS-TO for calling Snps and Indels + initial filter + split"
  }
  parameter_meta {
    input_bam: "Aligned Bam"
  }
  command <<<
    bash ~{m2_script} \
      --out-prefix "out/~{sample_name}~{suffix}" \
      --ref-fasta "~{ref_fasta}" \
      --input-bam "~{input_bam}" \
      --caller "~{mt_caller}" \
      --clair3-model "~{clair3_model}" \
      --clairs-platform "~{clairs_platform}" \
      --threads ~{n_cpu} \
      ~{"--interval-list " + mt_interval_list} \
      ~{"--min-af " + vaf_filter_threshold} \
      ~{"--blacklisted-sites " + blacklisted_sites} \
      ~{"--extra-args '" + m2_extra_args + "'"}
  >>>
  output {
    File raw_vcf = "out/~{sample_name}~{suffix}.raw.vcf"
    File raw_vcf_idx = "out/~{sample_name}~{suffix}.raw.vcf.idx"
    File stats = "out/~{sample_name}~{suffix}.raw.vcf.stats"
    File output_bamOut = "out/~{sample_name}~{suffix}.bamout.bam"

    File filtered_vcf = "out/~{sample_name}~{suffix}.filtered.vcf"
    File filtered_vcf_idx = "out/~{sample_name}~{suffix}.filtered.vcf.idx"

    File vcf_for_haplochecker = "out/~{sample_name}~{suffix}.splitAndPassOnly.vcf"
  }
}


###################################
# MongoM2FilterContaminationSplit #
####################################

task MongoM2FilterContaminationSplit {
  input {
    Int mem
    File raw_vcf
    File raw_vcf_stats
    String sample_name
    String hasContamination
    Float contamination_major
    Float contamination_minor
    Float nuc_contamination_estimate
    File filter_script
    File ref_fasta
    File ref_fai
    File ref_dict
    String suffix

    String? m2_extra_filtering_args
    Int max_alt_allele_count
    Float? vaf_filter_threshold
    Float? f_score_beta

    File? blacklisted_sites
    File? blacklisted_sites_index
  }
  Int command_mem = (mem * 1000) - 500

  meta {
    description: "Mutect2 Filtering for calling Snps and Indels"
  }
  parameter_meta {
      vaf_filter_threshold: "Hard cutoff for minimum allele fraction. All sites with VAF less than this cutoff will be filtered."
      f_score_beta: "F-Score beta balances the filtering strategy between recall and precision. The relative weight of recall to precision."
  }
  command <<<
    bash ~{filter_script} \
      --sample-name "~{sample_name}" \
      --suffix "~{suffix}" \
      --command-mem ~{command_mem} \
      --ref-fasta "~{ref_fasta}" \
      --raw-vcf "~{raw_vcf}" \
      --raw-vcf-stats "~{raw_vcf_stats}" \
      --has-contamination "~{hasContamination}" \
      --contamination-major ~{contamination_major} \
      --contamination-minor ~{contamination_minor} \
      --nuc-contamination-estimate ~{nuc_contamination_estimate} \
      --max-alt-allele-count ~{max_alt_allele_count} \
      ~{"--m2-extra-filtering-args '" + m2_extra_filtering_args + "'"} \
      ~{"--vaf-filter-threshold " + vaf_filter_threshold} \
      ~{"--f-score-beta " + f_score_beta} \
      ~{"--blacklisted-sites " + blacklisted_sites}
  >>>
  output {
    File filtered_vcf = "out/~{sample_name}~{suffix}.vcf"
    File filtered_vcf_idx = "out/~{sample_name}~{suffix}.vcf.idx"
    File split_vcf = "out/~{sample_name}~{suffix}.split.vcf"
    File split_vcf_index = "out/~{sample_name}~{suffix}.split.vcf.idx"
    Float contamination = read_float("out/~{sample_name}~{suffix}.hc_contam.txt") # now an optional output, producing UNDEFINED if not computed
  }
}


#########################
# MongoCallMtAndShifted #
#########################

task MongoCallMtAndShifted {
  input {
    Int mean_coverage
    String sample_base_name
    File mt_self
    File mt_self_index
    File mt_self_dict
    File input_bam
    File input_bai
    File mt_interval_list
    File force_call_vcf

    File shifted_mt_self
    File shifted_mt_self_index
    File shifted_mt_self_dict
    File shifted_input_bam
    File shifted_input_bai
    File shifted_mt_interval_list
    File shifted_force_call_vcf
    File call_mt_script

    String? m2_extra_args
    String? shifted_m2_extra_args

    Int max_reads_per_alignment_start = 75
    Boolean make_bamout = false
    String suffix
    String mt_caller = "clairs_to"
    String clair3_model = "r1041_e82_400bps_sup_v520"
    String clairs_platform = "ont_r10_dorado_sup_5khz"
    Int n_cpu = 2
    Float? vaf_filter_threshold
    Int mem
  }

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int command_mem = (mem * 1000) - 500

  meta {
    description: "Clair3 / ClairS-TO for calling Snps and Indels; runs on both shifted and nonshifted"
  }
  parameter_meta {
    input_bam: "Aligned Bam"
  }
  command <<<
    set -e
    call_one() { # PREFIX REF BAM INTERVALS FORCE EXTRA
      bash ~{call_mt_script} \
        --out-prefix "$1" \
        --ref-fasta "$2" \
        --input-bam "$3" \
        --interval-list "$4" \
        --force-call-vcf "$5" \
        --extra-args "$6" \
        --caller "~{mt_caller}" \
        --clair3-model "~{clair3_model}" \
        --clairs-platform "~{clairs_platform}" \
        --threads ~{n_cpu} \
        ~{"--min-af " + vaf_filter_threshold}
    }
    call_one "out/~{sample_base_name}~{suffix}" "~{mt_self}" "~{input_bam}" "~{mt_interval_list}" "~{force_call_vcf}" "~{default="" m2_extra_args}"
    call_one "out/~{sample_base_name}~{suffix}.shifted" "~{shifted_mt_self}" "~{shifted_input_bam}" "~{shifted_mt_interval_list}" "~{shifted_force_call_vcf}" "~{default="" shifted_m2_extra_args}"
  >>>
  output {
    File raw_vcf = "out/~{sample_base_name}~{suffix}.raw.vcf"
    File raw_vcf_idx = "out/~{sample_base_name}~{suffix}.raw.vcf.idx"
    File stats = "out/~{sample_base_name}~{suffix}.raw.vcf.stats"
    File output_bamOut = "out/~{sample_base_name}~{suffix}.bamout.bam"

    File shifted_raw_vcf = "out/~{sample_base_name}~{suffix}.shifted.raw.vcf"
    File shifted_raw_vcf_idx = "out/~{sample_base_name}~{suffix}.shifted.raw.vcf.idx"
    File shifted_stats = "out/~{sample_base_name}~{suffix}.shifted.raw.vcf.stats"
    File shifted_output_bamOut = "out/~{sample_base_name}~{suffix}.shifted.bamout.bam"
  }
}


##############################################
# MongoLiftoverCombineMergeFilterContamSplit #
##############################################

task MongoLiftoverCombineMergeFilterContamSplit {
  input {
    Int mem
    String sample_base_name
    String suffix

    File mt_self
    File mt_self_index
    File mt_self_dict
    File shifted_vcf
    File shifted_vcf_idx
    File non_shifted_vcf
    File non_shifted_vcf_idx
    File shifted_stats
    File non_shifted_stats
    File shift_back_chain

    String hasContamination
    Float contamination_major
    Float contamination_minor
    Float nuc_contamination_estimate
    File blacklisted_sites
    File blacklisted_sites_index
    String? m2_extra_filtering_args
    Int max_alt_allele_count
    Float? vaf_filter_threshold
    Float? f_score_beta
    String? tier_thresholds

  }
  Int command_mem = (mem * 1000) - 500
  String d = "$" # a stupid trick to get ${} indexing in bash to work in Cromwell

  meta {
    description: "Lifts over shifted vcf of control region and combines it with the rest of the chrM calls, then applies the LowAF / blacklist filters and splits multi-allelics."
  }
  parameter_meta {
    shifted_vcf: "VCF of control region on shifted reference"
    non_shifted_vcf: "VCF of the rest of chrM on original reference"
    mt_self: "Original (not shifted) chrM reference"
    shift_back_chain: "Chain file to lift over from shifted reference to original chrM"
  }
  command<<<
    set -e

    mkdir out

    this_sample=out/"~{sample_base_name}~{suffix}"
    this_shifted_vcf="~{shifted_vcf}"
    this_nonshifted_vcf="~{non_shifted_vcf}"
    this_ref_fasta="~{mt_self}"
    this_shiftback_chain="~{shift_back_chain}"
    this_shifted_stats="~{shifted_stats}"
    this_nonshifted_stats="~{non_shifted_stats}"
    this_blacklisted="~{blacklisted_sites}"

    this_has_contam="~{hasContamination}"
    this_nuc_contam="~{nuc_contamination_estimate}"
    this_contam_major="~{contamination_major}"
    this_contam_minor="~{contamination_minor}"

    if [[ "~{d}{this_has_contam}" == 'YES' ]]; then
      if (( $(echo "~{d}{this_contam_major} == 0.0"|bc -l) )); then
        this_hc_contamination="~{d}{this_contam_minor}"
      else
        this_hc_contamination=$( bc <<< "1-~{d}{this_contam_major}" )
      fi
    else
      this_hc_contamination=0.0
    fi

    echo "~{d}{this_hc_contamination}" > "~{d}{this_sample}.r2.hc_contam.txt"

    if (( $(echo "~{d}{this_nuc_contam} > ~{d}{this_hc_contamination}"|bc -l) )); then
      this_max_contamination="~{d}{this_nuc_contam}"
    else
      this_max_contamination="~{d}{this_hc_contamination}"
    fi

    echo "Nuclear estimate: ~{d}{this_nuc_contam}; HC: ~{d}{this_hc_contamination}; Max: ~{d}{this_max_contamination}"

    this_basename_raw="~{d}{this_sample}.raw"
    touch "~{d}{this_sample}.merge.bamout.bam"

    gatk LiftoverVcf \
      -I "~{d}{this_shifted_vcf}" \
      -O "~{d}{this_basename_raw}.shifted_back.vcf" \
      -R "~{d}{this_ref_fasta}" \
      -CHAIN "~{d}{this_shiftback_chain}" \
      -REJECT "~{d}{this_basename_raw}.rejected.vcf"

    cat "~{d}{this_nonshifted_vcf}" | perl -ne 'print unless /^##contig=<ID=(?!chrM)/' > chrMFiltVCF.vcf

    gatk MergeVcfs \
      -I "~{d}{this_basename_raw}.shifted_back.vcf" \
      -I chrMFiltVCF.vcf \
      -O "~{d}{this_basename_raw}.merged.vcf"

    cat "~{d}{this_basename_raw}.rejected.vcf" | grep ^chrM | wc -l > "~{d}{this_sample}.n_failed_vars.txt"

    # Clair has no Mutect2-style stats: keep the output contract with a placeholder.
    cat "~{d}{this_shifted_stats}" "~{d}{this_nonshifted_stats}" > "~{d}{this_basename_raw}.combined.stats"

    echo "Now filtering (LowAF, blacklist) and splitting..."
    cp "~{d}{this_basename_raw}.merged.vcf" step.vcf
    ~{if defined(vaf_filter_threshold) then 'bcftools filter -s LowAF -e "MAX(FMT/AF) < ' + select_first([vaf_filter_threshold]) + '" -Ov -o step_next.vcf step.vcf && mv step_next.vcf step.vcf' else ''}
    bcftools filter -s blacklisted_site -m + --mask-file "~{d}{this_blacklisted}" -Ov -o "~{d}{this_sample}.vcf" step.vcf
    gatk IndexFeatureFile -I "~{d}{this_sample}.vcf"

    bcftools norm -m -any -f "~{d}{this_ref_fasta}" -Ov -o split_all.vcf "~{d}{this_sample}.vcf"
    ~{if defined(vaf_filter_threshold) then 'bcftools filter -s LowAF -m + -e "FMT/AF < ' + select_first([vaf_filter_threshold]) + '" -Ov -o "' + d + '{this_sample}.split.vcf" split_all.vcf' else 'cp split_all.vcf "' + d + '{this_sample}.split.vcf"'}
    gatk IndexFeatureFile -I "~{d}{this_sample}.split.vcf"

    # Optional tiered heteroplasmy thresholds, e.g. 0.03,0.05,0.10 (PASS records only).
    ~{if defined(tier_thresholds) then 'IFS="," read -ra tiers <<< "' + select_first([tier_thresholds]) + '"' else 'tiers=()'}
    for t in "~{d}{tiers[@]}"; do
      tag=$(awk -v t="~{d}t" 'BEGIN {printf "%02d", t * 100 + 0.5}')
      bcftools view -f PASS,. -i "FMT/AF >= ~{d}t" -Ov -o "~{d}{this_sample}.tier~{d}{tag}.vcf" "~{d}{this_sample}.split.vcf"
    done
  >>>
  output {
    File rejected_raw_vcf = "out/~{sample_base_name}~{suffix}.raw.rejected.vcf"
    File merged_raw_vcf = "out/~{sample_base_name}~{suffix}.raw.merged.vcf"
    File merged_raw_vcf_idx = "out/~{sample_base_name}~{suffix}.raw.merged.vcf.idx"
    Int number_failed = read_int("out/~{sample_base_name}~{suffix}.n_failed_vars.txt")
    File raw_stats = "out/~{sample_base_name}~{suffix}.raw.combined.stats"

    File filtered_vcf = "out/~{sample_base_name}~{suffix}.vcf"
    File filtered_vcf_idx = "out/~{sample_base_name}~{suffix}.vcf.idx"
    File split_vcf = "out/~{sample_base_name}~{suffix}.split.vcf"
    File split_vcf_idx = "out/~{sample_base_name}~{suffix}.split.vcf.idx"
    Array[File] tier_vcfs = glob("out/~{sample_base_name}~{suffix}.tier*.vcf")
    Float contamination = read_float('out/~{sample_base_name}~{suffix}.r2.hc_contam.txt') # now an optional output, producing UNDEFINED if not computed
  }
}


#############################
# MongoProduceSelfReference #
#############################

task MongoProduceSelfReference {
  # updated 220616 to output a warning for overlapping mtDNA variants and just output the count, along
  # with a filtering pipeline to remove duplicates and prevent downstream errors
  input {
    Int hail_mem
    Int mem
    String sample_name
    File input_nuc_vcf
    File input_mt_vcf

    String suffix
    File ref_fasta
    File ref_fasta_index
    File nuc_interval_list
    File mt_ref_fasta
    File mt_ref_fasta_index
    File mt_interval_list
    File non_control_region_interval_list

    File fa_renaming_script
    File variant_bounds_script
    File check_hom_overlap_script
    File hail_force_call_script
    File self_reference_helper_script
    File produce_self_reference_script

    Int n_shift
  }
  Int command_mem = (mem * 1000) - 500

  String mt_intervals_basename = basename(mt_interval_list, ".interval_list")
  String noncontrol_basename = basename(non_control_region_interval_list, ".interval_list")
  String selfref_intervals_suffix = ".SelfRefLiftover.interval_list"

  command <<<
    bash ~{produce_self_reference_script} \
      --sample-name "~{sample_name}" \
      --suffix "~{suffix}" \
      --n-shift ~{n_shift} \
      --command-mem ~{command_mem} \
      --hail-mem ~{hail_mem} \
      --input-nuc-vcf "~{input_nuc_vcf}" \
      --input-mt-vcf "~{input_mt_vcf}" \
      --ref-fasta "~{ref_fasta}" \
      --nuc-interval-list "~{nuc_interval_list}" \
      --mt-ref-fasta "~{mt_ref_fasta}" \
      --mt-ref-fasta-index "~{mt_ref_fasta_index}" \
      --mt-interval-list "~{mt_interval_list}" \
      --non-control-region-interval-list "~{non_control_region_interval_list}" \
      --fa-renaming-script "~{fa_renaming_script}" \
      --variant-bounds-script "~{variant_bounds_script}" \
      --check-hom-overlap-script "~{check_hom_overlap_script}" \
      --hail-force-call-script "~{hail_force_call_script}" \
      --self-reference-helper-script "~{self_reference_helper_script}"
  >>>

  output {
    File self_fasta = "out/~{sample_name}~{suffix}.fasta"
    File self_fasta_index = "out/~{sample_name}~{suffix}.fasta.fai"
    File self_dict = "out/~{sample_name}~{suffix}.dict"

    File grch38_to_self_chain = "out/reference_to_~{sample_name}.chain"
    File filtered_vcf_ref_coord = "out/~{sample_name}.mtref.homoplasmies.vcf" # in reference coordinates

    File self_cat_fasta = "out/~{sample_name}andNuc~{suffix}.fasta"
    File self_cat_fasta_index = "out/~{sample_name}andNuc~{suffix}.fasta.fai"
    File self_cat_dict = "out/~{sample_name}andNuc~{suffix}.dict"

    File self_reference_nuc_fasta = "out/~{sample_name}NucOnly~{suffix}.fasta"
    File self_reference_nuc_fasta_index = "out/~{sample_name}NucOnly~{suffix}.fasta.fai"
    File self_reference_nuc_fasta_dict = "out/~{sample_name}NucOnly~{suffix}.dict"
    File grch38_to_self_nuc_chain = "out/reference_to_~{sample_name}NucOnly.chain"

    Int nuc_variants_dropped = read_int("out/~{sample_name}.nuc.removed.txt")
    Int mtdna_consensus_overlaps = read_int("out/~{sample_name}.mtdna_consensus_overlaps.txt")
    Int nuc_consensus_overlaps = read_int("out/~{sample_name}.nucdna_consensus_overlaps.txt")

    File self_shifted_fasta = "out/~{sample_name}~{suffix}.shifted_by_~{n_shift}_bases.fasta"
    File self_shifted_fasta_index = "out/~{sample_name}~{suffix}.shifted_by_~{n_shift}_bases.fasta.fai"
    File self_shifted_dict = "out/~{sample_name}~{suffix}.shifted_by_~{n_shift}_bases.dict"

    File shift_forward_chain = "out/~{sample_name}~{suffix}.shifted_by_~{n_shift}_bases.shift_fwd_~{n_shift}_bases.chain"
    File shift_back_chain = "out/~{sample_name}~{suffix}.shifted_by_~{n_shift}_bases.shift_back_~{n_shift}_bases.chain"

    File self_shifted_cat_fasta = "out/~{sample_name}~{suffix}.shifted_by_~{n_shift}_bases.cat.fasta"
    File self_shifted_cat_fasta_index = "out/~{sample_name}~{suffix}.shifted_by_~{n_shift}_bases.cat.fasta.fai"
    File self_shifted_cat_dict = "out/~{sample_name}~{suffix}.shifted_by_~{n_shift}_bases.cat.dict"

    File lifted_mt_intervals = "out/~{mt_intervals_basename}.~{sample_name}~{selfref_intervals_suffix}"
    File lifted_noncontrol_intervals = "out/~{noncontrol_basename}.~{sample_name}~{selfref_intervals_suffix}"
    File lifted_control_intervals = "out/control_region_shifted.chrM.~{sample_name}~{selfref_intervals_suffix}"

    File reversed_hom_vcf = "out/~{sample_name}~{suffix}.reversed.selfRef.homoplasmies.vcf.bgz"
    File reversed_hom_vcf_idx = "out/~{sample_name}~{suffix}.reversed.selfRef.homoplasmies.vcf.bgz.tbi"
    File reversed_hom_filters_vcf = "out/~{sample_name}~{suffix}.reversed.withfilters.selfRef.homoplasmies.vcf.bgz"
    File reversed_hom_filters_vcf_idx = "out/~{sample_name}~{suffix}.reversed.withfilters.selfRef.homoplasmies.vcf.bgz.tbi"
    File reversed_hom_vcf_shifted = "out/~{sample_name}~{suffix}.reversed.selfRef.shifted.homoplasmies.vcf.bgz"
    File reversed_hom_vcf_shifted_idx = "out/~{sample_name}~{suffix}.reversed.selfRef.shifted.homoplasmies.vcf.bgz.tbi"
  }

}


######################################
# MongoAlignToMtRegShiftedAndMetrics #
######################################

task MongoAlignToMtRegShiftedAndMetrics {
  input {
    Int mem
    File input_bam
    File coverage_metrics_script
    File align_script
    String aligner_commandline = "minimap2 -ax map-ont -Y --secondary=no -y $bash_ref_fasta"
    String sample_base_name
    String suffix

    File mt
    File mt_index
    File mt_dict

    File mt_cat
    File mt_cat_index
    File mt_cat_dict

    File mt_shifted
    File mt_shifted_index
    File mt_shifted_dict

    File mt_shifted_cat
    File mt_shifted_cat_index
    File mt_shifted_cat_dict

    File mt_interval_list

    String? read_name_regex

    Int? read_length
    Int? coverage_cap
    Int? n_cpu
    Int min_aligned_length = 0
    Float min_aligned_fraction = 0
  }
  Int command_mem = (mem * 1000) - 500

  Int this_cpu = select_first([n_cpu, 2])
  String this_aligner_commandline = aligner_commandline + " -t " + this_cpu

  Int read_length_for_optimization = select_first([read_length, 151])
  String this_output_bam_basename = basename(input_bam, '.bam') + ".remap" + suffix

  meta {
    description: "Aligns long reads with minimap2, filters NUMT-like reads and outputs coordinate sorted bams."
  }
  parameter_meta {
    input_bam: "Unmapped bam"
  }
  command <<<
    bash ~{align_script} \
      --sample-base-name "~{sample_base_name}" \
      --suffix "~{suffix}" \
      --command-mem ~{command_mem} \
      --read-length ~{read_length_for_optimization} \
      --aligner-commandline '~{this_aligner_commandline}' \
      --threads ~{this_cpu} \
      --min-aligned-length ~{min_aligned_length} \
      --min-aligned-fraction ~{min_aligned_fraction} \
      --input-bam "~{input_bam}" \
      --mt-cat "~{mt_cat}" \
      --mt-shifted-cat "~{mt_shifted_cat}" \
      --mt "~{mt}" \
      --mt-dict "~{mt_dict}" \
      --mt-shifted-dict "~{mt_shifted_dict}" \
      --mt-interval-list "~{mt_interval_list}" \
      --coverage-metrics-script "~{coverage_metrics_script}" \
      ~{"--coverage-cap " + coverage_cap}
  >>>
  output {
    File mt_aligned_bam = "out/~{this_output_bam_basename}.bam"
    File mt_aligned_bai = "out/~{this_output_bam_basename}.bai"
    File nuc_and_mt_aligned_bam = "out/~{this_output_bam_basename}_pre_mt_filt.bam"
    File nuc_and_mt_aligned_bai = "out/~{this_output_bam_basename}_pre_mt_filt.bai"
    File aligner_stderr_log = "out/~{this_output_bam_basename}.aligner.stderr.log"
    File duplicate_metrics = "out/~{this_output_bam_basename}.metrics"

    File wgs_metrics = "out/~{sample_base_name}~{suffix}_r2_wgs_metrics.txt"
    File theoretical_sensitivity = "out/~{sample_base_name}~{suffix}_r2_wgs_theoretical_sensitivity.txt"
    Int mean_coverage = read_int("out/~{sample_base_name}~{suffix}_r2_mean_coverage.txt")
    Float median_coverage = read_float("out/~{sample_base_name}~{suffix}_r2_median_coverage.txt")

    File shifted_mt_aligned_bam = "out/~{this_output_bam_basename}.shifted.bam"
    File shifted_mt_aligned_bai = "out/~{this_output_bam_basename}.shifted.bai"
    File nuc_and_shifted_mt_aligned_bam = "out/~{this_output_bam_basename}.shifted_pre_mt_filt.bam"
    File nuc_and_shifted_mt_aligned_bai = "out/~{this_output_bam_basename}.shifted_pre_mt_filt.bai"
    File shifted_aligner_stderr_log = "out/~{this_output_bam_basename}.shifted.aligner.stderr.log"
    File shifted_duplicate_metrics = "out/~{this_output_bam_basename}.shifted.metrics"
  }
}


##################################
# MongoLiftoverVCFAndGetCoverage #
##################################

task MongoLiftoverVCFAndGetCoverage {
  # A specialized routine to return the resultant VCF back to GRCh38
  input {
    Int hail_mem
    String sample_name
    File original_filtered_vcf
    File new_self_ref_vcf
    File reversed_hom_ref_vcf

    File mt_self
    File mt_self_index
    File mt_self_dict

    File mt_self_shifted
    File mt_self_shifted_index
    File mt_self_shifted_dict

    File chain_self_to_ref
    File chain_ref_to_self

    File input_bam_regular_ref
    File input_bam_regular_ref_index
    File input_bam_shifted_ref
    File input_bam_shifted_ref_index
    File self_control_region_shifted_reference_interval_list
    File self_non_control_region_interval_list

    String self_suffix
    File HailLiftover
    File liftover_vcf_script
    File coverage_helpers_script
    File ref_fasta
    File ref_fasta_index
    File ref_dict
  }
  command <<<
    bash ~{liftover_vcf_script} \
      --sample-name "~{sample_name}" \
      --self-suffix "~{self_suffix}" \
      --hail-mem ~{hail_mem} \
      --original-filtered-vcf "~{original_filtered_vcf}" \
      --new-self-ref-vcf "~{new_self_ref_vcf}" \
      --reversed-hom-ref-vcf "~{reversed_hom_ref_vcf}" \
      --mt-self "~{mt_self}" \
      --mt-self-index "~{mt_self_index}" \
      --mt-self-shifted "~{mt_self_shifted}" \
      --chain-self-to-ref "~{chain_self_to_ref}" \
      --chain-ref-to-self "~{chain_ref_to_self}" \
      --input-bam-regular-ref "~{input_bam_regular_ref}" \
      --input-bam-shifted-ref "~{input_bam_shifted_ref}" \
      --self-control-region-shifted-reference-interval-list "~{self_control_region_shifted_reference_interval_list}" \
      --self-non-control-region-interval-list "~{self_non_control_region_interval_list}" \
      --hail-liftover "~{HailLiftover}" \
      --coverage-helpers-script "~{coverage_helpers_script}" \
      --ref-fasta "~{ref_fasta}" \
      --ref-fasta-index "~{ref_fasta_index}"
  >>>


  output {
    File liftover_r1_rejected_vcf = "out/~{sample_name}~{self_suffix}.split.selfToRef.rejected.vcf"
    File liftover_r1_vcf = "out/~{sample_name}~{self_suffix}.split.selfToRef.pre.vcf"
    File liftover_r2_success_r1_vcf = "out/~{sample_name}~{self_suffix}.split.round2liftover.updated_success.vcf"
    File liftover_r2_rejected_vcf = "out/~{sample_name}~{self_suffix}.split.round2liftover.rejected.vcf"
    File liftover_r2_intermediate_vcf = "out/~{sample_name}~{self_suffix}.split.round2liftover.fixed.vcf"
    File liftover_r2_final_vcf = "out/~{sample_name}~{self_suffix}.split.selfToRef.final.vcf"
    File liftover_r2_log = "out/~{sample_name}~{self_suffix}.split_fix_liftover.log"
    File gap_coverage = "out/~{sample_name}~{self_suffix}.split.round2liftover.deletions_coverage.tsv"
    File self_coverage_table = "out/~{sample_name}~{self_suffix}.split.per_base_coverage.tsv"
    File liftoverStats = "out/~{sample_name}~{self_suffix}.split.round2liftover.all_int_outputs.final.txt"

    # stats
    Int n_liftover_changed_selfref_and_passed = read_int('out/~{sample_name}_n_ref_pass_thru.txt')
    Int n_liftover_r2_failed = read_int('out/~{sample_name}~{self_suffix}.split.round2liftover.round2_failed_sites.txt')
    Int n_liftover_r2_fixed = read_int('out/~{sample_name}~{self_suffix}.split.round2liftover.round2_fixed_sites.txt')
    Int n_liftover_r2_pass = read_int('out/~{sample_name}_n_final_pass.txt')
    Int n_liftover_r2_left_shift = read_int('out/~{sample_name}~{self_suffix}.split.round2liftover.left_alignment_of_indels.txt')
    Int n_liftover_r2_injected_from_success = read_int('out/~{sample_name}~{self_suffix}.split.round2liftover.round2_success_injected.txt')
    Int n_liftover_r2_ref_insertion_new_haplo = read_int('out/~{sample_name}~{self_suffix}.split.round2liftover.ref_insertion_new_haplos.txt')
    Int n_liftover_r2_failed_het_dele_span_insertion_boundary = read_int('out/~{sample_name}~{self_suffix}.split.round2liftover.het_deletions_span_insertions.txt')
    Int n_liftover_r2_failed_new_dupes_leftshift = read_int('out/~{sample_name}~{self_suffix}.split.round2liftover.new_dupes_left_shift_failed.txt')
    Int n_liftover_r2_het_ins_sharing_lhs_hom_dele = read_int('out/~{sample_name}~{self_suffix}.split.round2liftover.het_insertions_sharing_lhs_with_hom_ref_deletion.txt')
    Int n_liftover_r2_spanning_complex = read_int('out/~{sample_name}~{self_suffix}.split.round2liftover.het_dele_span_insert_repaired_with_complex_rework.txt')
    Int n_liftover_r2_spanningfixrhs_sharedlhs = read_int('out/~{sample_name}~{self_suffix}.split.round2liftover.heteroplasmic_deletions_sharing_lhs_with_homoplasmic_insertion_spanning_rhs.txt')
    Int n_liftover_r2_spanningfixlhs_upstream = read_int('out/~{sample_name}~{self_suffix}.split.round2liftover.heteroplasmic_deletions_spanning_only_lhs_of_homoplasmic_insertion.txt')
    Int n_liftover_r2_repaired_success = read_int('out/~{sample_name}~{self_suffix}.split.round2liftover.success_sites_flipped.txt')
  }
}


###############################
# CreateSpanIntervalsWithDict #
###############################

task CreateSpanIntervalsWithDict {
  input {
    File input_intervals
    File target_dict
    File self_reference_helper_script
    String intertext

  }
  String output_intervals = basename(input_intervals, ".interval_list") + intertext + ".SelfRefLiftover.interval_list"

  command <<<
    Rscript --vanilla ~{self_reference_helper_script} create-span-intervals \
      "~{input_intervals}" "~{target_dict}" "~{output_intervals}"
  >>>

  output {
    File lifted_intervals = "~{output_intervals}"
  }

}


###################
# MoveChainToZero #
###################

task MoveChainToZero {
  input {
    File source_chain
    File ref_intervals
    File self_reference_helper_script

  }
  String output_chain = basename(source_chain, ".chain") + ".toZero.chain"

  command <<<
    Rscript --vanilla ~{self_reference_helper_script} move-chain-to-zero \
      "~{source_chain}" "~{ref_intervals}" "~{output_chain}"
  >>>

  output {
    File chain = "~{output_chain}"
  }

}


##########################
# NucCoverageAtEveryBase #
##########################

task NucCoverageAtEveryBase {
  input {
    File input_bam_regular_ref
    File input_bam_regular_ref_index
    File input_bam_self_ref
    File input_bam_self_ref_index
    File input_bam_self_ref_shifted
    File input_bam_self_ref_shifted_index
    File nuc_interval_list
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File? self_nuc_interval_list
    File self_fasta
    File self_fasta_index
    File self_dict
    File? self_nuc_interval_list_shifted
    File self_shifted_fasta
    File self_shifted_fasta_index
    File self_shifted_dict

  }

  meta {
    description: "Mainly for QC to understand how NUMT coverage changes with remapping."
  }
  command <<<
    set -e

    gatk CollectHsMetrics \
      -I "~{input_bam_regular_ref}" \
      -R "~{ref_fasta}" \
      --PER_BASE_COVERAGE reference_pre_realignment.tsv \
      -O reference_pre_realignment.metrics \
      -TI "~{nuc_interval_list}" \
      -BI "~{nuc_interval_list}" \
      --COVERAGE_CAP 20000 \
      --SAMPLE_SIZE 1

    gatk CollectHsMetrics \
      -I "~{input_bam_self_ref}" \
      -R "~{self_fasta}" \
      --PER_BASE_COVERAGE self_post_realignment.tsv \
      -O self_post_realignment.metrics \
      -TI "~{self_nuc_interval_list}" \
      -BI "~{self_nuc_interval_list}" \
      --COVERAGE_CAP 20000 \
      --SAMPLE_SIZE 1

    gatk CollectHsMetrics \
      -I "~{input_bam_self_ref_shifted}" \
      -R "~{self_shifted_fasta}" \
      --PER_BASE_COVERAGE self_post_realignment_shifted.tsv \
      -O self_post_realignment_shifted.metrics \
      -TI "~{self_nuc_interval_list_shifted}" \
      -BI "~{self_nuc_interval_list_shifted}" \
      --COVERAGE_CAP 20000 \
      --SAMPLE_SIZE 1
  >>>
  output {
    File table_old_ref = "reference_pre_realignment.tsv"
    File table_new_self = "self_post_realignment.tsv"
    File table_new_self_shifted = "self_post_realignment_shifted.tsv"
  }
}


###########################
# MergeMultiSampleOutputs #
###########################


task MergeMultiSampleOutputs {
  # this task merges arrayed inputs for:
  # (1) output stats
  # (2) lifted-over variant calls
  # (3) per-base coverage
  input {
    Int hail_mem
    Array[String] sample_name
    Array[File] variant_vcf
    Array[File] coverage_table
    Array[File] statistics

    File MergePerBatch
    File coverage_helpers_script
  }

  command <<<
    set -e

    # start by merging statistics
    Rscript --vanilla ~{coverage_helpers_script} merge-statistics \
      "~{write_lines(statistics)}" batch_analysis_statistics.tsv

    # now produce the relevant inputs for the per-batch MT script
    Rscript --vanilla ~{coverage_helpers_script} create-batch-inputs \
      "~{write_lines(sample_name)}" \
      "~{write_lines(coverage_table)}" \
      "~{write_lines(variant_vcf)}" \
      coverage_paths.tsv \
      vcf_paths.tsv

    mkdir tmp

    # now merge coverage
    export PYSPARK_SUBMIT_ARGS="--driver-memory ~{hail_mem}g pyspark-shell"
    python ~{MergePerBatch} \
    --run-coverage \
    --input-tsv coverage_paths.tsv \
    --temp-dir tmp/ \
    --output-flat-file batch_merged_mt_coverage.tsv.bgz

    # finally merge variant calls
    python ~{MergePerBatch} \
    --run-variants \
    --input-tsv vcf_paths.tsv \
    --temp-dir tmp/ \
    --output-flat-file batch_merged_mt_calls.vcf.bgz
  >>>

  output {
    File merged_statistics = "batch_analysis_statistics.tsv"
    File merged_coverage = "batch_merged_mt_coverage.tsv.bgz"
    File merged_calls = "batch_merged_mt_calls.vcf.bgz"
  }

}


#############################
# MongoChainSwapLiftoverBed #
#############################

task MongoChainSwapLiftoverBed {
  input {
    File source_chain
    String input_target_name

    String input_source_name
    File input_bed
  }
  String this_bed_basename = basename(input_bed, '.bed')
  String d = "$" # a stupid trick to get ${} indexing in bash to work in Cromwell

  command <<<
    set -e

    mkdir out

    this_chain="~{source_chain}"
    this_target_sample="~{input_target_name}"
    this_target_path=out/"~{input_target_name}"
    this_base_path=out/"~{this_bed_basename}"

    chainSwap "~{d}{this_chain}" "~{d}{this_target_path}_to_~{input_source_name}.chain"
    liftOver ~{input_bed} "~{d}{this_chain}" "~{d}{this_base_path}.~{d}{this_target_sample}.liftedOver.bed" "~{d}{this_base_path}.~{d}{this_target_sample}.liftedOverRejects.unmapped"
    igvtools index "~{d}{this_base_path}.~{d}{this_target_sample}.liftedOver.bed"
  >>>

  output {
    File chain = "out/~{input_target_name}_to_~{input_source_name}.chain"
    File transformed_bed = "out/~{this_bed_basename}.~{input_target_name}.liftedOver.bed"
    File transformed_bed_index = "out/~{this_bed_basename}.~{input_target_name}.liftedOver.bed.idx"
    File rejected = "out/~{this_bed_basename}.~{input_target_name}.liftedOverRejects.unmapped"
  }

}


######################################
# MongoLiftoverSelfAndCollectOutputs #
######################################

task MongoLiftoverSelfAndCollectOutputs {
  # run liftover provided a coverage txt file
  # use reference allele depth from the SELF-REFERENCE data after reversing the force-call
  # to impute the appropriate coverage for a deleted site. This REF allele is the
  # actual original sequence.
  input {
    String sample_name
    File self_ref_table
    File chain
    File homoplasmic_deletions_coverage
    File liftover_table
    File coverage_helpers_script
    File collect_script

    Int mean_coverage
    Float median_coverage
    String major_haplogroup
    Float contamination
    Int nuc_variants_pass
    Int n_reads_unpaired_dropped
    Int nuc_variants_dropped
    Int mtdna_consensus_overlaps
    Int nuc_consensus_overlaps

  }

  command <<<
    bash ~{collect_script} \
      --sample-name "~{sample_name}" \
      --self-ref-table "~{self_ref_table}" \
      --chain "~{chain}" \
      --homoplasmic-deletions-coverage "~{homoplasmic_deletions_coverage}" \
      --liftover-table "~{liftover_table}" \
      --coverage-helpers-script "~{coverage_helpers_script}" \
      --mean-coverage ~{mean_coverage} \
      --median-coverage ~{median_coverage} \
      --major-haplogroup "~{major_haplogroup}" \
      --contamination ~{contamination} \
      --nuc-variants-pass ~{nuc_variants_pass} \
      --n-reads-unpaired-dropped ~{n_reads_unpaired_dropped} \
      --nuc-variants-dropped ~{nuc_variants_dropped} \
      --mtdna-consensus-overlaps ~{mtdna_consensus_overlaps} \
      --nuc-consensus-overlaps ~{nuc_consensus_overlaps}
  >>>

  output {
    File reference_coverage = "out/~{sample_name}.appended.liftedOver.tsv"
    File intermediate_bed = "out/~{sample_name}.liftedOver.bed"
    File rejected = "out/~{sample_name}.liftedOver.unmapped"
    File table = "out/~{sample_name}_mtanalysis_diagnostic_statistics.tsv"
  }

}


#############
# ChainSwap #
#############

task ChainSwap {
  # a shortcut method (that may be better behaved) which simply swaps reference and target
  input {
    File source_chain
    String input_source_name
    String input_target_name
  }
  String chain_output = input_target_name + "_to_" + input_source_name + ".chain"

  command {
    chainSwap ~{source_chain} ~{chain_output}
  }

  output {
    File chain = "~{chain_output}"
  }

}


###########################
# LiftOverAndJoinCoverage #
###########################

task LiftOverAndJoinCoverage {
  # As this is meant for QC, currently this does not allow spiked in
  # coverage at sites with gaps
  input {
    File ref_table
    File self_table
    File self_table_shifted
    File? chain
    File coverage_helpers_script
    Boolean? use_fast_join = true
  }

  String this_basename = basename(self_table, '.tsv')
  String output_bed = this_basename + ".nuc.liftedOver.bed"
  String output_bed_shifted = this_basename + ".shifted.nuc.liftedOver.bed"
  String rejects = this_basename + ".nuc.liftedOver.unmapped"
  String rejects_shifted = this_basename + ".shifted.nuc.liftedOver.unmapped"
  String final_tsv = this_basename + ".nuc.appended.liftedOver.tsv"

  command <<<
    tail -n +2 ~{self_table} | awk '{OFS="\t"} {print $1"\t"$2-1"\t"$2"\t"$4"\t"$3}' > per_base_coverage.bed
    tail -n +2 ~{self_table_shifted} | awk '{OFS="\t"} {print $1"\t"$2-1"\t"$2"\t"$4"\t"$3}' > per_base_coverage_shifted.bed
    tail -n +2 ~{ref_table} | awk '{OFS="\t"} {print $1"\t"$2-1"\t"$2"\t"$4"\t"$3}' > per_base_coverage_pre.bed
    liftOver per_base_coverage.bed ~{chain} ~{output_bed} ~{rejects}
    liftOver per_base_coverage_shifted.bed ~{chain} ~{output_bed_shifted} ~{rejects_shifted}

    Rscript --vanilla ~{coverage_helpers_script} join-lifted-coverage \
      "~{ref_table}" \
      per_base_coverage_pre.bed \
      "~{output_bed}" \
      "~{output_bed_shifted}" \
      cov.tsv \
      "~{use_fast_join}"

    bgzip -c cov.tsv > ~{final_tsv}.bgz
  >>>

  output {
    File reference_coverage = "~{final_tsv}.bgz"
    File intermediate_bed = "~{output_bed}"
    File rejected = "~{rejects}"
    File intermediate_bed_shifted = "~{output_bed_shifted}"
    File rejected_shifted = "~{rejects_shifted}"
  }

}


####################
# GetContamination #
####################

task GetContamination {
  input {
    Int mem
    File input_vcf
    File haplocheck_cli
    File contamination_script
    String sample_name
    Int mean_coverage
  }
  Int command_mem = (mem * 1000) - 500

  meta {
    description: "Uses new Haplochecker to estimate levels of contamination in mitochondria"
  }
  parameter_meta {
    input_vcf: "Filtered and split multi-allelic sites VCF for mitochondria"
  }
  command <<<
    bash ~{contamination_script} \
      --sample-name "~{sample_name}" \
      --command-mem ~{command_mem} \
      --mean-coverage ~{mean_coverage} \
      --input-vcf "~{input_vcf}" \
      --haplocheck-cli "~{haplocheck_cli}"
  >>>
  output {
    File contamination_file = "out/~{sample_name}_output_noquotes"
    String hasContamination = read_string("out/~{sample_name}.contamination.txt")
    String major_hg = read_string("out/~{sample_name}.major_hg.txt")
    String minor_hg = read_string("out/~{sample_name}.minor_hg.txt")
    Float major_level = read_float("out/~{sample_name}.mean_het_major.txt")
    Float minor_level = read_float("out/~{sample_name}.mean_het_minor.txt")
  }
}
