#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop(paste(
    "Expected an operation: stitch-self-coverage, merge-statistics,",
    "create-batch-inputs, append-lifted-coverage, write-diagnostic-statistics,",
    "or join-lifted-coverage."
  ))
}

operation <- args[1]
operation_args <- args[-1]

if (operation == "stitch-self-coverage") {
  if (length(operation_args) != 4) {
    stop("Usage: coverage_helpers.R stitch-self-coverage SELF_FASTA CONTROL_TSV NONCONTROL_TSV OUTPUT_TSV")
  }
  full_fasta <- readLines(operation_args[1])
  nlen <- nchar(paste0(full_fasta[2:length(full_fasta)], collapse = ""))
  nshift <- 8000
  shift_back <- function(x) {
    if (x < (nlen - nshift + 1)) {
      return(x + nshift)
    }
    x - (nlen - nshift)
  }

  control_region_shifted <- read.table(operation_args[2], header = TRUE)
  control_region_shifted[, "pos"] <- sapply(control_region_shifted[, "pos"], shift_back)
  beginning <- subset(control_region_shifted, control_region_shifted[, "pos"] < 8000)
  end <- subset(control_region_shifted, control_region_shifted[, "pos"] > 8000)
  non_control_region <- read.table(operation_args[3], header = TRUE)
  combined_table <- rbind(beginning, non_control_region, end)
  write.table(
    combined_table,
    operation_args[4],
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    sep = "\t"
  )
} else if (operation == "merge-statistics") {
  if (length(operation_args) != 2) {
    stop("Usage: coverage_helpers.R merge-statistics STATISTICS_PATHS OUTPUT_TSV")
  }
  files_of_interest <- read.csv(
    operation_args[1],
    sep = "\t",
    stringsAsFactors = FALSE,
    header = FALSE
  )[[1]]
  dfs <- lapply(files_of_interest, function(path) {
    read.csv(path, sep = "\t", stringsAsFactors = FALSE)
  })
  df <- do.call("rbind", dfs)
  write.table(
    df,
    sep = "\t",
    row.names = FALSE,
    file = operation_args[2],
    quote = FALSE
  )
} else if (operation == "create-batch-inputs") {
  if (length(operation_args) != 5) {
    stop("Usage: coverage_helpers.R create-batch-inputs SAMPLES COVERAGE_PATHS VCF_PATHS COVERAGE_OUTPUT VCF_OUTPUT")
  }
  sample_ids <- read.csv(
    operation_args[1],
    sep = "\t",
    stringsAsFactors = FALSE,
    header = FALSE
  )[[1]]
  paths_coverage <- read.csv(
    operation_args[2],
    sep = "\t",
    stringsAsFactors = FALSE,
    header = FALSE
  )[[1]]
  paths_vcf <- read.csv(
    operation_args[3],
    sep = "\t",
    stringsAsFactors = FALSE,
    header = FALSE
  )[[1]]
  write.table(
    data.frame(s = sample_ids, path = paths_coverage),
    sep = "\t",
    row.names = FALSE,
    file = operation_args[4],
    quote = FALSE
  )
  write.table(
    data.frame(s = sample_ids, path = paths_vcf),
    sep = "\t",
    row.names = FALSE,
    file = operation_args[5],
    quote = FALSE
  )
} else if (operation == "append-lifted-coverage") {
  if (length(operation_args) != 5) {
    stop("Usage: coverage_helpers.R append-lifted-coverage CHAIN FAILED_BED SUCCESS_BED DELETION_COVERAGE OUTPUT_TSV")
  }
  init_chain <- readLines(operation_args[1])
  if (length(grep("^chain", init_chain)) != 1) {
    stop("ERROR: only single block chain files are currently supported.")
  }
  loc_init_chain <- grep("^chain", init_chain)
  chain_splits <- strsplit(init_chain[loc_init_chain], "\\s+")
  if (!all(sapply(chain_splits, function(x) as.numeric(c(x[6], x[11]))) == 0)) {
    stop("ERROR: chain file must start at 0 for both source and target.")
  }
  if (!(1 %in% loc_init_chain)) {
    stop("ERROR: first line of chain file must be header.")
  }
  len_tab <- length(strsplit(init_chain[2], "\t")[[1]])
  len_spc <- length(strsplit(init_chain[2], " ")[[1]])
  sep <- if (len_tab > len_spc) "\t" else " "

  chain <- read.csv(
    operation_args[1],
    skip = 1,
    header = FALSE,
    col.names = c("bases_in_block", "ns", "nt"),
    sep = sep,
    stringsAsFactors = FALSE,
    fill = TRUE
  )
  chain[nrow(chain), 2:3] <- 0
  if (nrow(chain) == 1) {
    chain[["coord_s"]] <- chain[["bases_in_block"]]
    chain[["coord_t"]] <- chain[["bases_in_block"]]
  } else {
    chain[["coord_s"]] <- cumsum(chain[["bases_in_block"]]) +
      cumsum(c(0, chain[["ns"]][2:nrow(chain) - 1]))
    chain[["coord_t"]] <- cumsum(chain[["bases_in_block"]]) +
      cumsum(c(0, chain[["nt"]][2:nrow(chain) - 1]))
  }
  if (any(is.na(chain))) {
    stop("ERROR: there should be no missing fields in chain.")
  }
  if (sum((chain[["ns"]] > 0) & (chain[["nt"]] > 0)) != 0) {
    stop(paste(
      "ERROR: Chain files should not have breaks -- e.g., each ungapped block",
      "should be followed by a gap in EITHER the source or target, but not both."
    ))
  }

  start_end_to_bed <- function(start, end_inclusive) {
    new_start <- start + 1
    ends <- end_inclusive
    if (length(new_start) != length(ends)) {
      stop("ERROR: starts and ends must have the same length.")
    }
    if (length(new_start) == 0) {
      empty_df <- data.frame("chr" = "chrM", "start" = 1, "end" = 1)
      return(empty_df[empty_df[["start"]] > 1, ])
    }
    bp <- unlist(lapply(seq_along(start), function(idx) new_start[idx]:ends[idx]))
    df <- data.frame("chr" = "chrM", "start" = bp - 1, "end" = bp)
    df[order(df[["start"]]), ]
  }

  bed_entries_missing_in_db <- function(query, db) {
    if ((nrow(query) > 0) & (nrow(db) > 0)) {
      query[
        !paste0(query[["start"]], ":", query[["end"]]) %in%
          paste0(db[["start"]], ":", db[["end"]]),
      ]
    } else {
      query
    }
  }
  n_rows_different_keys_bed <- function(df1, df2) {
    nrow(bed_entries_missing_in_db(df1, df2)) +
      nrow(bed_entries_missing_in_db(df2, df1))
  }

  insertions <- chain[chain[["ns"]] > 0, ]
  expected_failed_liftover <- start_end_to_bed(
    insertions[["coord_s"]],
    insertions[["coord_s"]] + insertions[["ns"]]
  )
  failure_bed <- read.csv(
    operation_args[2],
    header = FALSE,
    sep = "\t",
    col.names = c("chr", "start", "end", "cov"),
    comment.char = "#"
  )
  if (n_rows_different_keys_bed(expected_failed_liftover, failure_bed) > 0) {
    stop("ERROR: Rejected bed file is not identical to expectation.")
  }

  target_len <- as.numeric(strsplit(init_chain[1], "\\s+")[[1]][9])
  deletions <- chain[chain[["nt"]] > 0, ]
  expected_liftover_gaps <- start_end_to_bed(
    deletions[["coord_t"]],
    deletions[["coord_t"]] + deletions[["nt"]]
  )
  success_bed <- read.csv(
    operation_args[3],
    header = FALSE,
    sep = "\t",
    col.names = c("chr", "start", "end", "cov"),
    comment.char = "#"
  )
  full_bed <- data.frame("chr" = "chrM", "start" = 0:(target_len - 1), "end" = 1:target_len)
  not_in_success_bed <- bed_entries_missing_in_db(full_bed, success_bed)
  if (n_rows_different_keys_bed(expected_liftover_gaps, not_in_success_bed) > 0) {
    stop("ERROR: Success bed file has different gaps than expected based on chain.")
  }
  if (nrow(expected_liftover_gaps) == 0) {
    print("NOTE: No gaps found in coverage file, so no changes were made.")
    appended_bed <- success_bed[order(success_bed[["start"]]), ]
  } else {
    fill_in_coverage <- read.csv(operation_args[4], sep = "\t", stringsAsFactors = FALSE)
    bed_homoplas_dele <- data.frame(
      "chr" = "chrM",
      "start" = fill_in_coverage[["reference_position"]] - 1,
      "end" = fill_in_coverage[["reference_position"]],
      "cov" = fill_in_coverage[["ref_allele_depth"]]
    )
    if (n_rows_different_keys_bed(expected_liftover_gaps, bed_homoplas_dele) > 0) {
      stop("ERROR: Records in homoplasmic_deletions_coverage are different from expected based on chain.")
    }
    appended_bed <- rbind(
      success_bed,
      merge(not_in_success_bed, bed_homoplas_dele, by = c("chr", "start", "end"))
    )
    appended_bed <- appended_bed[order(appended_bed[["start"]]), ]
    row.names(appended_bed) <- NULL
    if ((nrow(bed_entries_missing_in_db(full_bed, appended_bed)) > 0) |
        (!all(appended_bed[["end"]] == 1:target_len))) {
      stop("ERROR: the final lifted over bed file should have one row per position on the mtDNA.")
    }
  }
  final_tsv <- appended_bed
  names(final_tsv) <- c("chrom", "to_rm", "pos", "coverage")
  final_tsv[["target"]] <- "."
  write.table(
    final_tsv[, c("chrom", "pos", "target", "coverage")],
    operation_args[5],
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )
} else if (operation == "write-diagnostic-statistics") {
  if (length(operation_args) != 11) {
    stop("Usage: coverage_helpers.R write-diagnostic-statistics LIFTOVER_STATS OUTPUT MEAN MEDIAN HAPLOGROUP CONTAMINATION NUC_PASS READS_DROPPED NUC_DROPPED MT_OVERLAPS NUC_OVERLAPS")
  }
  liftover_stats <- read.csv(operation_args[1], sep = "\t", stringsAsFactors = FALSE)
  metrics <- data.frame(
    mean_coverage = operation_args[3],
    median_coverage = operation_args[4],
    major_haplogroup = operation_args[5],
    contamination = operation_args[6],
    nuc_variants_pass = operation_args[7],
    n_reads_unpaired_dropped = operation_args[8],
    nuc_variants_dropped = operation_args[9],
    mtdna_consensus_overlaps = operation_args[10],
    nuc_consensus_overlaps = operation_args[11]
  )
  write.table(
    cbind(liftover_stats, metrics),
    operation_args[2],
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )
} else if (operation == "join-lifted-coverage") {
  if (length(operation_args) != 6) {
    stop("Usage: coverage_helpers.R join-lifted-coverage REF_TABLE REF_BED SELF_BED SHIFTED_BED OUTPUT USE_FAST_JOIN")
  }
  use_fast_join <- operation_args[6] == "true"
  ref_pre_bed <- read.csv(
    operation_args[2],
    header = FALSE,
    sep = "\t",
    col.names = c("chr", "start", "end", "cov_ref", "target"),
    comment.char = "#",
    stringsAsFactors = FALSE
  )
  self_lifted_bed <- read.csv(
    operation_args[3],
    header = FALSE,
    sep = "\t",
    col.names = c("chr", "start", "end", "cov_self_ref", "target"),
    comment.char = "#",
    stringsAsFactors = FALSE
  )
  self_lifted_shifted_bed <- read.csv(
    operation_args[4],
    header = FALSE,
    sep = "\t",
    col.names = c("chr", "start", "end", "cov_self_ref_shifted", "target"),
    comment.char = "#",
    stringsAsFactors = FALSE
  )
  original_cov <- read.csv(operation_args[1], sep = "\t", stringsAsFactors = FALSE)

  edit_self_bed <- function(self_lifted_bed) {
    chr_to_convert <- unique(self_lifted_bed[["chr"]])
    mapped_chrs <- lapply(chr_to_convert, function(x) {
      unique(original_cov[["chrom"]][original_cov[["target"]] == x])
    })
    if (!all(sapply(mapped_chrs, length) == 1)) {
      stop(paste(
        "ERROR: some intervals did not map to chromosomes in the original",
        "reference coverage files (or mapped to multiple)."
      ))
    }
    chr_mapper <- sapply(mapped_chrs, function(x) x[[1]])
    names(chr_mapper) <- chr_to_convert
    self_lifted_bed[["chr"]] <- chr_mapper[self_lifted_bed[["chr"]]]
    split_by_targ <- split(self_lifted_bed, self_lifted_bed[["target"]])
    split_by_targ_updated <- lapply(names(split_by_targ), function(x) {
      this_tab <- split_by_targ[[x]]
      filtered_ref_start <- ref_pre_bed[ref_pre_bed[["target"]] == x, ][["start"]]
      this_tab[["start"]] <- this_tab[["start"]] + filtered_ref_start[1]
      this_tab[["end"]] <- this_tab[["end"]] + filtered_ref_start[1]
      this_tab
    })
    do.call("rbind", split_by_targ_updated)
  }

  ensure_unique <- function(tab, keep_cols = FALSE) {
    tf_redundant <- all(tab[["end"]] - tab[["start"]]) == 1
    if (tf_redundant) {
      tab[["joiner"]] <- paste0(tab[["chr"]], ":", tab[["start"]], ":", tab[["target"]])
      vals <- c("chr", "start", "target")
    } else {
      stop("ERROR: end is not the same as start+1 -- this is not handled.")
      tab[["joiner"]] <- paste0(
        tab[["chr"]], ":", tab[["start"]], ":", tab[["end"]], ":",
        self_lifted_bed[["target"]]
      )
      vals <- c("chr", "start", "end", "target")
    }
    if (!keep_cols) {
      tab <- tab[, -1 * which(names(tab) %in% vals)]
    }
    if (any(duplicated(tab[["joiner"]]))) {
      stop("Joining rows not unique. This means you have duplicate loci in the per-base coverage results.")
    }
    list("vals" = vals, "tab" = tab)
  }

  fast_join <- function(tab1, tab2, by, target) {
    all_keys <- unique(c(tab1[[by]], tab2[[by]]))
    new_keys <- all_keys[!all_keys %in% tab1[[by]]]
    if (length(new_keys) > 0) {
      new_df <- data.frame(ff = new_keys)
      names(new_df) <- by
      new_df[setdiff(names(tab1), names(new_df))] <- NA
      tab1 <- rbind(tab1, new_df)
    }
    tab1[[target]] <- tab2[[target]][match(tab1[[by]], tab2[[by]])]
    tab1
  }

  self_lifted_bed <- edit_self_bed(self_lifted_bed)
  self_lifted_shifted_bed <- edit_self_bed(self_lifted_shifted_bed)
  if (use_fast_join) {
    ref_pre_bed <- ensure_unique(ref_pre_bed, TRUE)
    self_lifted_bed <- ensure_unique(self_lifted_bed)
    self_lifted_shifted_bed <- ensure_unique(self_lifted_shifted_bed)
    joint_bed <- fast_join(
      ref_pre_bed[["tab"]],
      self_lifted_bed[["tab"]],
      by = "joiner",
      target = "cov_self_ref"
    )
    joint_bed <- fast_join(
      joint_bed,
      self_lifted_shifted_bed[["tab"]],
      by = "joiner",
      target = "cov_self_ref_shifted"
    )
    names(joint_bed) <- c(
      "chrom", "to_rm", "pos", "coverage_original", "target", "to_rm2",
      "coverage_remapped_self", "coverage_remapped_self_shifted"
    )
  } else {
    joint_bed <- merge(
      ref_pre_bed,
      self_lifted_bed,
      by = c("chr", "start", "end", "target"),
      all = TRUE
    )
    joint_bed <- merge(
      joint_bed,
      self_lifted_shifted_bed,
      by = c("chr", "start", "end", "target"),
      all = TRUE
    )
    names(joint_bed) <- c(
      "chrom", "to_rm", "pos", "target", "coverage_original",
      "coverage_remapped_self", "coverage_remapped_self_shifted"
    )
  }
  write.table(
    joint_bed[, c(
      "chrom", "pos", "target", "coverage_original",
      "coverage_remapped_self", "coverage_remapped_self_shifted"
    )],
    operation_args[5],
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )
} else {
  stop(paste0("Unknown operation: ", operation))
}
