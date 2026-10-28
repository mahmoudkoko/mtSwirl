#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop(paste(
    "Expected an operation: shift-reference, adjust-mt-intervals,",
    "shift-control-region, create-span-intervals, or move-chain-to-zero."
  ))
}

operation <- args[1]
operation_args <- args[-1]

if (operation == "shift-reference") {
  if (length(operation_args) != 5) {
    stop("Usage: self_reference_helpers.R shift-reference FASTA N_SHIFT SHIFTED_FASTA BACK_CHAIN FORWARD_CHAIN")
  }
  full_fasta <- readLines(operation_args[1])
  topline <- full_fasta[1]
  linelen <- nchar(full_fasta[2])
  n_shift <- as.numeric(operation_args[2])
  other_lines <- paste0(full_fasta[2:length(full_fasta)], collapse = "")
  other_lines_shifted <- paste0(
    substr(other_lines, n_shift + 1, nchar(other_lines)),
    substr(other_lines, 1, n_shift)
  )
  len_chr <- nchar(other_lines_shifted)
  shifted_split <- substring(
    other_lines_shifted,
    seq(1, len_chr, linelen),
    unique(c(seq(linelen, len_chr, linelen), len_chr))
  )
  writeLines(c(topline, shifted_split), operation_args[3])

  total_len <- nchar(other_lines)
  sec1 <- paste(
    c("chain", 9999, "chrM", total_len, "+", 0, total_len - n_shift,
      "chrM", total_len, "+", n_shift, total_len, 1),
    collapse = " "
  )
  sec2 <- paste(
    c("chain", 9999, "chrM", total_len, "+", total_len - n_shift, total_len,
      "chrM", total_len, "+", 0, n_shift, 2),
    collapse = " "
  )
  writeLines(c(sec1, total_len - n_shift, "", sec2, n_shift, ""), operation_args[4])

  sec1_f <- paste(
    c("chain", 9999, "chrM", total_len, "+", n_shift, total_len,
      "chrM", total_len, "+", 0, total_len - n_shift, 1),
    collapse = " "
  )
  sec2_f <- paste(
    c("chain", 9999, "chrM", total_len, "+", 0, n_shift,
      "chrM", total_len, "+", total_len - n_shift, total_len, 2),
    collapse = " "
  )
  writeLines(c(sec1_f, total_len - n_shift, "", sec2_f, n_shift, ""), operation_args[5])
} else if (operation == "adjust-mt-intervals") {
  if (length(operation_args) != 3) {
    stop("Usage: self_reference_helpers.R adjust-mt-intervals INTERVAL_LIST DICT OUTPUT")
  }
  intervals <- readLines(operation_args[1])
  intervals <- intervals[grep("^@", intervals, invert = TRUE)]
  interval_names <- sapply(strsplit(intervals, "\t"), function(x) x[5])
  new_header <- readLines(operation_args[2])
  lens <- sapply(interval_names, function(x) {
    dict_row <- new_header[grep(paste0("^@SQ\tSN:", x), new_header)[1]]
    as.numeric(gsub("LN:", "", strsplit(dict_row, "\t")[[1]][3]))
  })
  if (any(is.na(lens))) {
    stop("ERROR: Some NUMT intervals were not found in the mt_andNuc sequence dictionary.")
  }
  new_intervals <- c(
    new_header,
    paste(interval_names, 1, lens, "+", interval_names, sep = "\t")
  )
  writeLines(new_intervals, operation_args[3])
} else if (operation == "shift-control-region") {
  if (length(operation_args) != 4) {
    stop("Usage: self_reference_helpers.R shift-control-region NONCONTROL_INTERVAL_LIST SHIFTED_DICT N_SHIFT OUTPUT")
  }
  full_intervals <- readLines(operation_args[1])
  correct_dict <- readLines(operation_args[2])
  n_shift <- as.numeric(operation_args[3])
  if (length(full_intervals) > 3) {
    stop("ERROR: there should be only 3 lines in interval list.")
  }
  if ((length(grep("^@", full_intervals)) != 2) ||
      (length(grep("^chrM", full_intervals)) != 1)) {
    stop("ERROR: there should be 2 comment lines and one interval line on chrM")
  }
  out_intervals <- c(full_intervals[1], correct_dict[2])
  split_line2 <- strsplit(out_intervals[2], "\t")[[1]]
  this_len <- as.numeric(gsub("^LN:", "", split_line2[grep("^LN", split_line2)[1]]))

  interval <- strsplit(full_intervals[3], "\t")[[1]]
  new_end <- as.numeric(interval[2]) - 1 + this_len - n_shift
  new_start <- as.numeric(interval[3]) + 1 - n_shift
  interval[2:3] <- c(as.character(new_start), as.character(new_end))
  out_intervals <- c(out_intervals, paste0(interval, collapse = "\t"))
  writeLines(out_intervals, operation_args[4])
} else if (operation == "create-span-intervals") {
  if (length(operation_args) != 3) {
    stop("Usage: self_reference_helpers.R create-span-intervals INTERVALS DICT OUTPUT")
  }
  intervals <- readLines(operation_args[1])
  intervals <- intervals[grep("^@", intervals, invert = TRUE)]
  interval_names <- sapply(strsplit(intervals, "\t"), function(x) x[5])
  new_header <- readLines(operation_args[2])
  lens <- sapply(interval_names, function(x) {
    dict_row <- new_header[grep(paste0("^@SQ\tSN:", x), new_header)[1]]
    as.numeric(gsub("LN:", "", strsplit(dict_row, "\t")[[1]][3]))
  })
  if (any(is.na(lens))) {
    stop("ERROR: Some NUMT intervals were not found in the mt_andNuc sequence dictionary.")
  }
  new_intervals <- c(
    new_header,
    paste(interval_names, 1, lens, "+", interval_names, sep = "\t")
  )
  writeLines(new_intervals, operation_args[3])
} else if (operation == "move-chain-to-zero") {
  if (length(operation_args) != 3) {
    stop("Usage: self_reference_helpers.R move-chain-to-zero CHAIN INTERVALS OUTPUT")
  }
  chain <- readLines(operation_args[1])
  intervals <- readLines(operation_args[2])
  intervals <- intervals[grep("^@", intervals, invert = TRUE)]
  intervals_split <- strsplit(intervals, "\t")
  target_names <- sapply(intervals_split, function(x) x[5])
  names(target_names) <- sapply(
    intervals_split,
    function(x) paste0(x[1], ":", as.numeric(x[2]) - 1, "-", x[3])
  )
  to_edit <- grep("^chain", chain)
  new_headers <- sapply(chain[to_edit], function(x) {
    this_header <- strsplit(x, " ")[[1]]
    this_search <- paste0(this_header[8], ":", this_header[11], "-", this_header[12])
    if (!this_search %in% names(target_names)) {
      stop("ERROR: chain file has a segment not in intervals list.")
    }
    this_name <- target_names[this_search]
    this_header[c(3, 8)] <- this_name
    this_header[c(4, 6, 7)] <- as.numeric(this_header[c(4, 6, 7)]) -
      as.numeric(this_header[6])
    this_header[c(9, 11, 12)] <- as.numeric(this_header[c(9, 11, 12)]) -
      as.numeric(this_header[11])
    paste0(this_header, collapse = " ")
  })
  new_chain <- chain
  new_chain[to_edit] <- new_headers
  writeLines(new_chain, operation_args[3])
} else {
  stop(paste0("Unknown operation: ", operation))
}
