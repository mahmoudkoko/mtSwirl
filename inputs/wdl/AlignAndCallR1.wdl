version 1.0

import "MongoTasks.wdl" as MongoTasks_Single

workflow AlignAndCallR1 {
  meta {
    description: "Takes in unmapped bam and outputs VCF of SNP/Indel calls on the mitochondria."
  }

  input {
    File NucHaplotypeCallerScript
    File Mutect2NucScript
    File Mutect2MtInitialScript
    File FilterContaminationScript
    File GetContaminationScript
    File input_bam
    File input_bai
    String sample_name

    File ref_dict
    File ref_fasta
    File ref_fasta_index

    File mt_dict
    File mt_fasta
    File mt_fasta_index
    File blacklisted_sites
    File blacklisted_sites_index

    File nuc_interval_list
    File mt_interval_list

    Int mt_mean_coverage

    Boolean use_haplotype_caller_nucdna
    Int haplotype_caller_nucdna_dp_lower_bound
    String? m2_extra_args
    String? m2_filter_extra_args
    Float? vaf_filter_threshold
    Float? f_score_beta
    Float nuc_contamination_estimate

    File haplocheck_cli

    String mt_caller = "clairs_to"
    String clair3_model = "r1041_e82_400bps_sup_v520"
    String clairs_platform = "ont_r10_dorado_sup_5khz"
    Int n_cpu = 2

    Int mem_small
    Int mem_standard
    Int mem_mutect_mt
  }

  parameter_meta {
  }

  if (use_haplotype_caller_nucdna) {
    call MongoTasks_Single.MongoHC as CallNucHCIntegrated {
      input:
        hc_script = NucHaplotypeCallerScript,
        input_bam = input_bam,
        input_bai = input_bai,
        sample_name = sample_name,
        nuc_interval_list = nuc_interval_list,
        ref_fasta = ref_fasta,
        ref_fai = ref_fasta_index,
        ref_dict = ref_dict,
        suffix = '.nuc',
        haplotype_caller_nucdna_dp_lower_bound = haplotype_caller_nucdna_dp_lower_bound,
        mt_caller = mt_caller,
        clair3_model = clair3_model,
        clairs_platform = clairs_platform,
        n_cpu = n_cpu,
        mem = mem_standard
    }
  }
  if (!use_haplotype_caller_nucdna) {
    call MongoTasks_Single.MongoNucM2 as CallNucM2Integrated {
      input:
        m2_script = Mutect2NucScript,
        input_bam = input_bam,
        input_bai = input_bai,
        sample_name = sample_name,

        ref_fasta = ref_fasta,
        ref_fai = ref_fasta_index,
        ref_dict = ref_dict,
        suffix = '.nuc',
        mt_interval_list = nuc_interval_list,

        m2_extra_args = select_first([m2_extra_args, ""]),

        max_alt_allele_count = 4,
        vaf_filter_threshold = 0.95,
        mt_caller = mt_caller,
        clair3_model = clair3_model,
        clairs_platform = clairs_platform,
        n_cpu = n_cpu,
        nuc_contamination_estimate = nuc_contamination_estimate,

        mem = mem_standard
    }
  }

  call MongoTasks_Single.MongoRunM2InitialFilterSplit as CallMt {
    input:
      m2_script = Mutect2MtInitialScript,
      sample_name = sample_name,
      input_bam = input_bam,
      input_bai = input_bai,
      nuc_contamination_estimate = nuc_contamination_estimate,
      mt_interval_list = mt_interval_list,
      ref_fasta = mt_fasta,
      ref_fai = mt_fasta_index,
      ref_dict = mt_dict,
      suffix = "",
      m2_extra_filtering_args = select_first([m2_filter_extra_args, ""]) + " --min-median-mapping-quality 0",
      max_alt_allele_count = 4,
      vaf_filter_threshold = select_first([vaf_filter_threshold, 0.005]),
      mt_caller = mt_caller,
      clair3_model = clair3_model,
      clairs_platform = clairs_platform,
      n_cpu = n_cpu,
      blacklisted_sites = blacklisted_sites,
      blacklisted_sites_index = blacklisted_sites_index,
      f_score_beta = f_score_beta,
      m2_extra_args = select_first([m2_extra_args, ""]),
      mean_coverage = mt_mean_coverage,
      mem = mem_mutect_mt
  }

  call MongoTasks_Single.GetContamination {
    input:
      contamination_script = GetContaminationScript,
      input_vcf = CallMt.vcf_for_haplochecker,
      sample_name = sample_name,
      mean_coverage = mt_mean_coverage,
      haplocheck_cli = haplocheck_cli,
      mem = mem_small
  }

  call MongoTasks_Single.MongoM2FilterContaminationSplit as FilterContamination {
    input:
      filter_script = FilterContaminationScript,
      raw_vcf = CallMt.filtered_vcf,
      raw_vcf_stats = CallMt.stats,
      sample_name = sample_name,
      hasContamination = GetContamination.hasContamination,
      contamination_major = GetContamination.major_level,
      contamination_minor = GetContamination.minor_level,
      suffix = "",
      nuc_contamination_estimate = nuc_contamination_estimate,
      ref_fasta = mt_fasta,
      ref_fai = mt_fasta_index,
      ref_dict = mt_dict,
      m2_extra_filtering_args = select_first([m2_filter_extra_args, ""]) + " --min-median-mapping-quality 0",
      max_alt_allele_count = 4,
      vaf_filter_threshold = vaf_filter_threshold,
      blacklisted_sites = blacklisted_sites,
      blacklisted_sites_index = blacklisted_sites_index,
      f_score_beta = f_score_beta,
      mem = mem_standard
  }

  output {
    File out_vcf = FilterContamination.filtered_vcf
    File out_vcf_index = FilterContamination.filtered_vcf_idx
    File split_vcf = FilterContamination.split_vcf
    File split_vcf_index = FilterContamination.split_vcf_index
    File nuc_vcf = select_first([CallNucHCIntegrated.full_pass_vcf, CallNucM2Integrated.full_pass_vcf])
    File nuc_vcf_index = select_first([CallNucHCIntegrated.full_pass_vcf_index, CallNucM2Integrated.full_pass_vcf_index])
    File nuc_vcf_unfiltered = select_first([CallNucHCIntegrated.filtered_vcf, CallNucM2Integrated.filtered_vcf])
    File split_nuc_vcf = select_first([CallNucHCIntegrated.split_vcf, CallNucM2Integrated.split_vcf])
    File split_nuc_vcf_index = select_first([CallNucHCIntegrated.split_vcf_index, CallNucM2Integrated.split_vcf_index])
    Int nuc_variants_pass = select_first([CallNucHCIntegrated.post_filt_vars, CallNucM2Integrated.post_filt_vars])
    File input_vcf_for_haplochecker = CallMt.vcf_for_haplochecker
    File contamination_metrics = GetContamination.contamination_file
    String major_haplogroup = GetContamination.major_hg
    Float contamination = FilterContamination.contamination
    String hasContamination = GetContamination.hasContamination
    Float contamination_major = GetContamination.major_level
    Float contamination_minor = GetContamination.minor_level
  }
}
