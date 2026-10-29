version 1.0

import "ProcessSample.wdl" as ProcessSingleSample
import "MongoTasks.wdl" as MongoTasks

workflow SwirlPipeline {

  meta {
    description: "Takes in a list of hg38 bam or cram and outputs VCF of SNP/Indel calls on the mitochondria."
    allowNestedInputs: true
  }

  input {
    File wgs_aligned_input_bam_or_cram_list
    File wgs_aligned_input_bam_or_cram_index_list
    File sample_name_list

    File mt_interval_list
    File nuc_interval_list

    # Read length used for optimization only. If this is too small CollectWgsMetrics might fail, but the results are not
    # affected by this number. Default is 151.
    Int? max_read_length

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File mt_dict
    File mt_fasta
    File mt_fasta_index
    File blacklisted_sites
    File blacklisted_sites_index

    File non_control_region_interval_list

    File HailLiftover
    File LiftoverVcfScript
    File HailForceCallScript
    File CoverageMetricsScript
    File AlignToMtRefScript
    File SubsetBamScript
    File NucHaplotypeCallerScript
    File CallMtScript
    File Mutect2NucScript
    File Mutect2MtInitialScript
    File FilterContaminationScript
    File GetContaminationScript
    File LiftoverCollectScript
    File SelfReferenceHelperScript
    File ProduceSelfReferenceScript
    File CoverageHelpersScript
    File haplocheck_cli
    File FaRenamingScript
    File CheckVariantBoundsScript
    File CheckHomOverlapScript
    File MergePerBatch

    String? requester_pays_project
    String? m2_extra_args
    String? m2_filter_extra_args
    String? printreads_extra_args
    Float? vaf_filter_threshold
    Float? f_score_beta
    Float nuc_contamination_estimate = 0.0
    Boolean compute_numt_coverage = false
    Boolean use_haplotype_caller_nucdna = true
    Boolean skip_restore_hardclips = false
    Int haplotype_caller_nucdna_dp_lower_bound = 10
    # Memory in GB per group of tasks; each task's JVM heap is mem*1000 - 500 MB.
    Int mem_small = 2
    Int mem_standard = 4
    Int mem_mutect_mt = 6
    Int mem_align = 6
    Int hail_mem = 4
    Int? n_cpu_bwa
  }

  Array[File] input_cram = read_lines(wgs_aligned_input_bam_or_cram_list)
  Array[File] input_crai = read_lines(wgs_aligned_input_bam_or_cram_index_list)
  Array[String] sample_name = read_lines(sample_name_list)

  scatter (cram_pair in zip(zip(input_cram, input_crai), sample_name)) {
    call ProcessSingleSample.MitochondriaPipeline as ProcessSingleSample {
        input:
          wgs_aligned_input_bam_or_cram = cram_pair.left.left,
          wgs_aligned_input_bam_or_cram_index = cram_pair.left.right,
          sample_name = cram_pair.right,
          mt_interval_list = mt_interval_list,
          nuc_interval_list = nuc_interval_list,

          max_read_length = max_read_length,
          skip_restore_hardclips = skip_restore_hardclips,

          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_dict = ref_dict,

          mt_dict = mt_dict,
          mt_fasta = mt_fasta,
          mt_fasta_index = mt_fasta_index,
          blacklisted_sites = blacklisted_sites,
          blacklisted_sites_index = blacklisted_sites_index,

          non_control_region_interval_list = non_control_region_interval_list,

          HailLiftover = HailLiftover,
          LiftoverVcfScript = LiftoverVcfScript,
          HailForceCallScript = HailForceCallScript,
          SelfReferenceHelperScript = SelfReferenceHelperScript,
          ProduceSelfReferenceScript = ProduceSelfReferenceScript,
          CoverageHelpersScript = CoverageHelpersScript,
          CoverageMetricsScript = CoverageMetricsScript,
          AlignToMtRefScript = AlignToMtRefScript,
          SubsetBamScript = SubsetBamScript,
          NucHaplotypeCallerScript = NucHaplotypeCallerScript,
          CallMtScript = CallMtScript,
          Mutect2NucScript = Mutect2NucScript,
          Mutect2MtInitialScript = Mutect2MtInitialScript,
          FilterContaminationScript = FilterContaminationScript,
          GetContaminationScript = GetContaminationScript,
          LiftoverCollectScript = LiftoverCollectScript,
          haplocheck_cli = haplocheck_cli,
          FaRenamingScript = FaRenamingScript,
          CheckVariantBoundsScript = CheckVariantBoundsScript,
          CheckHomOverlapScript = CheckHomOverlapScript,

          requester_pays_project = requester_pays_project,
          m2_extra_args = m2_extra_args,
          m2_filter_extra_args = m2_filter_extra_args,
          printreads_extra_args = printreads_extra_args,
          vaf_filter_threshold = vaf_filter_threshold,
          f_score_beta = f_score_beta,
          nuc_contamination_estimate = nuc_contamination_estimate,
          compute_numt_coverage = compute_numt_coverage,
          use_haplotype_caller_nucdna = use_haplotype_caller_nucdna,
          haplotype_caller_nucdna_dp_lower_bound = haplotype_caller_nucdna_dp_lower_bound,

          mem_small = mem_small,
          mem_standard = mem_standard,
          mem_mutect_mt = mem_mutect_mt,
          mem_align = mem_align,
          hail_mem = hail_mem,
          n_cpu_bwa = n_cpu_bwa
    }
  }

  call MongoTasks.MergeMultiSampleOutputs {
    input:
      sample_name = sample_name,
      variant_vcf = ProcessSingleSample.final_vcf,
      coverage_table = ProcessSingleSample.final_base_level_coverage_metrics,
      statistics = ProcessSingleSample.stats_outputs,
      MergePerBatch = MergePerBatch,
      coverage_helpers_script = CoverageHelpersScript,
      hail_mem = hail_mem
  }

  output {
    Array[File] subset_bam = ProcessSingleSample.subset_bam
    Array[File] subset_bai = ProcessSingleSample.subset_bai
    Array[File] r1_vcf = ProcessSingleSample.r1_vcf
    Array[File] r1_vcf_index = ProcessSingleSample.r1_vcf_index
    Array[File] r1_nuc_vcf = ProcessSingleSample.r1_nuc_vcf
    Array[File] r1_nuc_vcf_index = ProcessSingleSample.r1_nuc_vcf_index
    Array[File] r1_nuc_vcf_unfiltered = ProcessSingleSample.r1_nuc_vcf_unfiltered
    Array[File] r1_split_vcf = ProcessSingleSample.r1_split_vcf
    Array[File] r1_split_vcf_index = ProcessSingleSample.r1_split_vcf_index

    Array[File] self_mt_aligned_bam = ProcessSingleSample.self_mt_aligned_bam
    Array[File] self_mt_aligned_bai = ProcessSingleSample.self_mt_aligned_bai
    Array[File] self_ref_vcf = ProcessSingleSample.self_ref_vcf
    Array[File] self_ref_vcf_index = ProcessSingleSample.self_ref_vcf_index
    Array[File] self_ref_split_vcf = ProcessSingleSample.self_ref_split_vcf
    Array[File] self_ref_split_vcf_index = ProcessSingleSample.self_ref_split_vcf_index
    Array[File] self_base_level_coverage_metrics = ProcessSingleSample.self_base_level_coverage_metrics
    Array[File] self_reference_fasta = ProcessSingleSample.self_reference_fasta
    Array[File] reference_to_self_ref_chain = ProcessSingleSample.reference_to_self_ref_chain
    Array[File] self_control_region_shifted = ProcessSingleSample.self_control_region_shifted
    Array[File] self_non_control_region = ProcessSingleSample.self_non_control_region

    Array[File] liftover_fix_pipeline_log = ProcessSingleSample.liftover_fix_pipeline_log
    Array[File] stats_outputs = ProcessSingleSample.stats_outputs

    Array[File] final_vcf = ProcessSingleSample.final_vcf
    Array[File] final_rejected_vcf = ProcessSingleSample.final_rejected_vcf
    Array[File] final_base_level_coverage_metrics = ProcessSingleSample.final_base_level_coverage_metrics
    Array[File?] numt_base_level_coverage = ProcessSingleSample.numt_base_level_coverage
    
    Array[File] input_vcf_for_haplochecker = ProcessSingleSample.input_vcf_for_haplochecker
    Array[File] duplicate_metrics = ProcessSingleSample.duplicate_metrics
    Array[File] coverage_metrics = ProcessSingleSample.coverage_metrics
    Array[File] theoretical_sensitivity_metrics = ProcessSingleSample.theoretical_sensitivity_metrics
    Array[File] contamination_metrics = ProcessSingleSample.contamination_metrics

    # liftover stats
    Array[Int] success_liftover_variants = ProcessSingleSample.success_liftover_variants
    Array[Int] failed_liftover_variants = ProcessSingleSample.failed_liftover_variants
    Array[Int] fixed_liftover_variants = ProcessSingleSample.fixed_liftover_variants
    Array[Int] n_liftover_r2_left_shift = ProcessSingleSample.n_liftover_r2_left_shift
    Array[Int] n_liftover_r2_injected_from_success = ProcessSingleSample.n_liftover_r2_injected_from_success
    Array[Int] n_liftover_r2_ref_insertion_new_haplo = ProcessSingleSample.n_liftover_r2_ref_insertion_new_haplo
    Array[Int] n_liftover_r2_failed_het_dele_span_insertion_boundary = ProcessSingleSample.n_liftover_r2_failed_het_dele_span_insertion_boundary
    Array[Int] n_liftover_r2_failed_new_dupes_leftshift = ProcessSingleSample.n_liftover_r2_failed_new_dupes_leftshift
    Array[Int] n_liftover_r2_het_ins_sharing_lhs_hom_dele = ProcessSingleSample.n_liftover_r2_het_ins_sharing_lhs_hom_dele
    Array[Int] n_liftover_r2_spanning_complex = ProcessSingleSample.n_liftover_r2_spanning_complex
    Array[Int] n_liftover_r2_spanningfixrhs_sharedlhs = ProcessSingleSample.n_liftover_r2_spanningfixrhs_sharedlhs
    Array[Int] n_liftover_r2_spanningfixlhs_upstream = ProcessSingleSample.n_liftover_r2_spanningfixlhs_upstream
    Array[Int] n_liftover_r2_repaired_success = ProcessSingleSample.n_liftover_r2_repaired_success

    # other stats
    Array[Int] mean_coverage = ProcessSingleSample.mean_coverage
    Array[Float] median_coverage = ProcessSingleSample.median_coverage
    Array[String] major_haplogroup = ProcessSingleSample.major_haplogroup
    Array[Float?] contamination = ProcessSingleSample.contamination
    Array[Int] nuc_variants_pass = ProcessSingleSample.nuc_variants_pass
    Array[Int] n_reads_unpaired_dropped = ProcessSingleSample.n_reads_unpaired_dropped
    Array[Int] nuc_variants_dropped = ProcessSingleSample.nuc_variants_dropped
    Array[Int] mtdna_consensus_overlaps = ProcessSingleSample.mtdna_consensus_overlaps
    Array[Int] nuc_consensus_overlaps = ProcessSingleSample.nuc_consensus_overlaps

    # merged files
    File merged_statistics = MergeMultiSampleOutputs.merged_statistics
    File merged_coverage = MergeMultiSampleOutputs.merged_coverage
    File merged_calls = MergeMultiSampleOutputs.merged_calls
  }
}
