#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: extract_wgs_coverage_metrics.R METRICS MEAN_OUTPUT MEDIAN_OUTPUT")
}

metrics <- read.table(
  args[1],
  skip = 6,
  header = TRUE,
  stringsAsFactors = FALSE,
  sep = "\t",
  nrows = 1
)
write.table(
  floor(metrics[, "MEAN_COVERAGE"]),
  args[2],
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE
)
write.table(
  metrics[, "MEDIAN_COVERAGE"],
  args[3],
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE
)
