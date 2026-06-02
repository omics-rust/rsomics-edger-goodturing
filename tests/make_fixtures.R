#!/usr/bin/env Rscript
# Generate a small golden count matrix (deterministic) for compat fixtures.
set.seed(42)
n_genes <- 400L
n_samp <- 4L
# Counts with lots of zeros/low values so the FoF table is rich (NB mixture).
m <- matrix(0L, n_genes, n_samp)
for (j in 1:n_samp) {
  lambda <- rgamma(n_genes, shape = 0.4, rate = 0.15)
  m[, j] <- as.integer(rpois(n_genes, lambda))
}
rownames(m) <- sprintf("gene%04d", 1:n_genes)
colnames(m) <- sprintf("s%d", 1:n_samp)
write.table(m, file = commandArgs(trailingOnly = TRUE)[1],
            sep = "\t", quote = FALSE, col.names = NA)
