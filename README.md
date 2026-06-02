# rsomics-edger-goodturing

Simple Good-Turing frequency smoothing of integer count vectors into proportion
estimates, including the mass reserved for unseen items. A Rust port of edgeR's
`goodTuring()` / `goodTuringProportions()`.

```bash
cargo install rsomics-edger-goodturing
```

## Usage

Input is a tab-separated integer count matrix (genes × samples), first column =
gene id, first row = sample names.

```bash
# per-distinct-count smoothed proportions + P0 per column
rsomics-edger-goodturing counts.tsv -o gt.tsv

# same-shape proportion matrix (goodTuringProportions): zeros share the unseen
# mass P0/n0, each column renormalized to sum 1
rsomics-edger-goodturing counts.tsv --proportions -o prop.tsv
```

`--conf` (default 1.96) sets the confidence multiplier for the Turing → linear
Good-Turing switch rule.

## What it computes

For each column, Simple Good-Turing builds the frequency-of-frequencies table
N_r (how many items were seen exactly r times), then:

- `P0 = N_1 / N` is the total probability reserved for unseen items.
- Each observed count r gets a smoothed proportion via the Turing estimate
  `r* = (r+1)·N_{r+1}/N_r`, switching to the log-linear `Z_r` regression
  estimate once the two agree within `conf·σ`.
- Proportions are renormalized so the seen mass sums to `1 − P0`.

A column with a single distinct count cannot be regressed and yields `NaN`,
matching edgeR.

## Origin

This crate is an independent Rust reimplementation of edgeR's `goodTuring`
based on:

- Gale, W. A. & Sampson, G. (1995). *Good-Turing frequency estimation without
  tears.* Journal of Quantitative Linguistics 2(3): 217–237.
  DOI: [10.1080/09296179508590051](https://doi.org/10.1080/09296179508590051)
- The public behaviour of `edgeR::goodTuring` / `edgeR::goodTuringProportions`
  observed black-box (output structure, zero handling, renormalization).

edgeR's Simple Good-Turing core is GPL-licensed C++; no source from it was used
during implementation. The algorithm here is reconstructed from the Gale &
Sampson paper and verified against the upstream binary's output. The edgeR R
wrappers (which only marshal the per-column data) informed the matrix-level
zero/renormalization conventions.

License: MIT OR Apache-2.0.
Upstream credit: edgeR <https://bioconductor.org/packages/edgeR/> (GPL ≥ 2);
McCarthy, Chen & Smyth (2012), Nucleic Acids Research 40(10): 4288–4297.
