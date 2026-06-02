#!/usr/bin/env Rscript
# edgeR goodTuring oracle. Args: <counts.tsv> <mode: prop|gt> <conf>
suppressMessages(library(edgeR))
a <- commandArgs(trailingOnly = TRUE)
counts <- as.matrix(read.table(a[1], header = TRUE, row.names = 1, sep = "\t", check.names = FALSE))
mode <- a[2]
conf <- as.numeric(a[3])

fmt <- function(v) ifelse(is.nan(v), "NaN", sprintf("%.10f", v))

if (mode == "prop") {
  z <- goodTuringProportions(counts)
  cat("gene", colnames(counts), sep = "\t"); cat("\n")
  for (i in seq_len(nrow(z))) {
    cat(rownames(counts)[i], fmt(z[i, ]), sep = "\t"); cat("\n")
  }
} else {
  cat("sample\tr\tNr\tproportion\tP0\tn0\n")
  for (j in seq_len(ncol(counts))) {
    g <- goodTuring(counts[, j], conf = conf)
    for (k in seq_along(g$count)) {
      cat(colnames(counts)[j], g$count[k], g$n[k], fmt(g$proportion[k]),
          fmt(g$P0), g$n0, sep = "\t"); cat("\n")
    }
  }
}
