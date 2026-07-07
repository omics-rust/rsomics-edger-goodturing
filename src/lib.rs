use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

pub struct Matrix {
    pub samples: Vec<String>,
    pub genes: Vec<String>,
    pub counts: Vec<i64>,
    pub n_samples: usize,
}

impl Matrix {
    pub fn load(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
        let mut lines = BufReader::new(file).lines();

        let header = lines
            .next()
            .ok_or_else(|| RsomicsError::InvalidInput("empty count matrix".into()))?
            .map_err(RsomicsError::Io)?;
        let samples: Vec<String> = header.split('\t').skip(1).map(str::to_string).collect();
        let n_samples = samples.len();

        let mut genes = Vec::new();
        let mut counts = Vec::new();
        for line in lines {
            let line = line.map_err(RsomicsError::Io)?;
            if line.is_empty() {
                continue;
            }
            let mut fields = line.split('\t');
            let gene = fields
                .next()
                .ok_or_else(|| RsomicsError::InvalidInput("row without a gene id".into()))?;
            genes.push(gene.to_string());
            let before = counts.len();
            for f in fields {
                let c = f.parse::<i64>().map_err(|_| {
                    RsomicsError::InvalidInput(format!("non-integer count '{f}' for gene {gene}"))
                })?;
                if c < 0 {
                    return Err(RsomicsError::InvalidInput(format!(
                        "negative count {c} for gene {gene}"
                    )));
                }
                counts.push(c);
            }
            if counts.len() - before != n_samples {
                return Err(RsomicsError::InvalidInput(format!(
                    "gene {gene}: {} values, header has {n_samples} samples",
                    counts.len() - before
                )));
            }
        }
        Ok(Self {
            samples,
            genes,
            counts,
            n_samples,
        })
    }
}

/// Result of Simple Good-Turing smoothing on one count vector.
pub struct GoodTuring {
    /// Distinct positive counts r, ascending.
    pub count: Vec<i64>,
    /// Frequency of frequency N_r at each r.
    pub n: Vec<i64>,
    /// Number of zero-count items.
    pub n0: i64,
    /// Smoothed proportion for an item observed r times; same length as `count`.
    pub proportion: Vec<f64>,
    /// Total mass reserved for unseen items, N_1 / N.
    pub p0: f64,
}

/// Simple Good-Turing smoothing (Gale & Sampson 1995). `x` is a count vector;
/// zeros are tallied into `n0` but not smoothed. Mirrors edgeR::goodTuring.
pub fn good_turing(x: &[i64], conf: f64) -> GoodTuring {
    let mut max_x = 0i64;
    for &v in x {
        if v > max_x {
            max_x = v;
        }
    }

    // edgeR picks the tabulation strategy by max(x) vs length(x): a dense
    // frequency table costs O(max_x), so when a single count dwarfs the vector
    // (e.g. a 1e10 library size) it would allocate gigabytes. In that regime we
    // instead collect frequency-of-frequencies from the sorted distinct counts,
    // which costs O(n log n) in the vector length. Both paths yield identical
    // (r, N_r) pairs, so the downstream math is unchanged.
    let (n0, count, nr) = if max_x < x.len() as i64 {
        let mut tab = vec![0i64; (max_x + 1) as usize];
        for &v in x {
            tab[v as usize] += 1;
        }
        let mut count = Vec::new();
        let mut nr = Vec::new();
        for (r, &c) in tab.iter().enumerate().skip(1) {
            if c > 0 {
                count.push(r as i64);
                nr.push(c);
            }
        }
        (tab[0], count, nr)
    } else {
        let mut sorted: Vec<i64> = x.to_vec();
        sorted.sort_unstable();
        let mut n0 = 0i64;
        let mut count = Vec::new();
        let mut nr = Vec::new();
        for &v in &sorted {
            if v == 0 {
                n0 += 1;
                continue;
            }
            if count.last() == Some(&v) {
                *nr.last_mut().unwrap() += 1;
            } else {
                count.push(v);
                nr.push(1);
            }
        }
        (n0, count, nr)
    };

    let proportion = simple_good_turing(&count, &nr, conf);
    let total: i64 = count.iter().zip(&nr).map(|(&r, &n)| r * n).sum();
    let p0 = if total == 0 || count[0] != 1 {
        0.0
    } else {
        nr[0] as f64 / total as f64
    };

    GoodTuring {
        count,
        n: nr,
        n0,
        proportion,
        p0,
    }
}

/// Core SGT: given distinct counts `r` and frequencies-of-frequencies `nr`,
/// return the smoothed proportion per r. A single distinct count cannot be
/// regressed, so it yields NaN — matching edgeR.
fn simple_good_turing(r: &[i64], nr: &[i64], conf: f64) -> Vec<f64> {
    let k = r.len();
    if k == 0 {
        return Vec::new();
    }
    if k == 1 {
        return vec![f64::NAN];
    }

    let total: i64 = r.iter().zip(nr).map(|(&ri, &ni)| ri * ni).sum();
    let p0 = if r[0] == 1 { nr[0] as f64 } else { 0.0 } / total as f64;

    // Z_r averages each N_r over the gap to its neighbours (paper eq. for Z).
    let mut log_r = Vec::with_capacity(k);
    let mut log_z = Vec::with_capacity(k);
    for j in 0..k {
        let i = if j == 0 { 0 } else { r[j - 1] };
        let kk = if j == k - 1 { 2 * r[j] - i } else { r[j + 1] };
        let z = 2.0 * nr[j] as f64 / (kk - i) as f64;
        log_r.push((r[j] as f64).ln());
        log_z.push(z.ln());
    }

    let (a, b) = linear_fit(&log_r, &log_z);
    let sr = |x: f64| (a + b * x.ln()).exp();

    let mut rstar = vec![0.0f64; k];
    let mut use_lgt = false;
    for j in 0..k {
        let rj = r[j] as f64;
        let lgt = (rj + 1.0) * sr(rj + 1.0) / sr(rj);
        if use_lgt {
            rstar[j] = lgt;
            continue;
        }
        // Turing estimate needs N_{r+1}; absent (or past the switch) → LGT.
        let next = if j + 1 < k && r[j + 1] == r[j] + 1 {
            Some(nr[j + 1] as f64)
        } else {
            None
        };
        match next {
            Some(np1) => {
                let nj = nr[j] as f64;
                let turing = (rj + 1.0) * np1 / nj;
                let sigma = ((rj + 1.0).powi(2) * (np1 / (nj * nj)) * (1.0 + np1 / nj)).sqrt();
                if (turing - lgt).abs() <= conf * sigma {
                    use_lgt = true;
                    rstar[j] = lgt;
                } else {
                    rstar[j] = turing;
                }
            }
            None => {
                use_lgt = true;
                rstar[j] = lgt;
            }
        }
    }

    let n_prime: f64 = nr.iter().zip(&rstar).map(|(&n, &rs)| n as f64 * rs).sum();
    rstar.iter().map(|&rs| (1.0 - p0) * rs / n_prime).collect()
}

fn linear_fit(x: &[f64], y: &[f64]) -> (f64, f64) {
    let n = x.len() as f64;
    let mean_x = x.iter().sum::<f64>() / n;
    let mean_y = y.iter().sum::<f64>() / n;
    let mut sxy = 0.0;
    let mut sxx = 0.0;
    for (&xi, &yi) in x.iter().zip(y) {
        sxy += (xi - mean_x) * (yi - mean_y);
        sxx += (xi - mean_x) * (xi - mean_x);
    }
    let slope = sxy / sxx;
    (mean_y - slope * mean_x, slope)
}

pub struct Opts {
    pub proportions: bool,
    pub conf: f64,
}

pub fn run(counts_path: &Path, opts: &Opts, output: &mut dyn Write) -> Result<usize> {
    let m = Matrix::load(counts_path)?;
    let n_genes = m.genes.len();
    if m.n_samples == 0 {
        return Err(RsomicsError::InvalidInput(
            "matrix has no sample columns".into(),
        ));
    }

    let mut out = BufWriter::new(output);
    if opts.proportions {
        write_proportions(&m, opts.conf, &mut out)?;
    } else {
        write_per_count(&m, opts.conf, &mut out)?;
    }
    out.flush().map_err(RsomicsError::Io)?;
    Ok(n_genes)
}

// Per-element proportion matrix, same shape as input: each entry's smoothed
// proportion, with zeros sharing the unseen mass P0/n0 — edgeR
// goodTuringProportions. Columns renormalize to sum 1.
fn write_proportions(m: &Matrix, conf: f64, out: &mut dyn Write) -> Result<()> {
    let n = m.genes.len();
    let s = m.n_samples;
    let mut prop = vec![0.0f64; n * s];

    let mut col = Vec::with_capacity(n);
    for j in 0..s {
        col.clear();
        col.extend((0..n).map(|i| m.counts[i * s + j]));
        let gt = good_turing(&col, conf);
        let per_zero = if gt.n0 > 0 { gt.p0 / gt.n0 as f64 } else { 0.0 };
        for i in 0..n {
            let c = col[i];
            prop[i * s + j] = if c == 0 {
                per_zero
            } else {
                let idx = gt.count.binary_search(&c).unwrap();
                gt.proportion[idx]
            };
        }
    }

    out.write_all(b"gene").map_err(RsomicsError::Io)?;
    for label in &m.samples {
        write!(out, "\t{label}").map_err(RsomicsError::Io)?;
    }
    out.write_all(b"\n").map_err(RsomicsError::Io)?;
    for (i, gene) in m.genes.iter().enumerate() {
        out.write_all(gene.as_bytes()).map_err(RsomicsError::Io)?;
        for j in 0..s {
            write!(out, "\t{}", fmt(prop[i * s + j])).map_err(RsomicsError::Io)?;
        }
        out.write_all(b"\n").map_err(RsomicsError::Io)?;
    }
    Ok(())
}

// Per-distinct-count table (goodTuring): r, N_r, smoothed proportion, plus a
// P0 summary line, per sample column.
fn write_per_count(m: &Matrix, conf: f64, out: &mut dyn Write) -> Result<()> {
    let n = m.genes.len();
    let s = m.n_samples;
    writeln!(out, "sample\tr\tNr\tproportion\tP0\tn0").map_err(RsomicsError::Io)?;
    let mut col = Vec::with_capacity(n);
    for j in 0..s {
        col.clear();
        col.extend((0..n).map(|i| m.counts[i * s + j]));
        let gt = good_turing(&col, conf);
        for ((&r, &nr), &p) in gt.count.iter().zip(&gt.n).zip(&gt.proportion) {
            writeln!(
                out,
                "{}\t{r}\t{nr}\t{}\t{}\t{}",
                m.samples[j],
                fmt(p),
                fmt(gt.p0),
                gt.n0
            )
            .map_err(RsomicsError::Io)?;
        }
    }
    Ok(())
}

fn fmt(v: f64) -> String {
    if v.is_nan() {
        "NaN".to_string()
    } else {
        format!("{v:.10}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn paper_example_matches_gale_sampson() {
        // Gale & Sampson 1995 Table-1 r/Nr, expanded back into a count vector.
        let r = [1i64, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 26];
        let nr = [120i64, 40, 24, 13, 15, 5, 11, 2, 2, 1, 3, 1];
        let mut x = Vec::new();
        for (&ri, &ni) in r.iter().zip(&nr) {
            for _ in 0..ni {
                x.push(ri);
            }
        }
        let gt = good_turing(&x, 1.96);
        assert_eq!(gt.count, r);
        assert!((gt.p0 - 0.196078431373).abs() < 1e-9);
        let expect = [
            0.001174424833,
            0.002042632783,
            0.003589686771,
            0.005220115910,
            0.006893291373,
            0.008591266915,
            0.010304901932,
            0.012029053987,
            0.013760609240,
            0.015497572281,
            0.018982778748,
            0.043536426607,
        ];
        for (got, want) in gt.proportion.iter().zip(&expect) {
            assert!((got - want).abs() < 1e-9, "got {got} want {want}");
        }
        let sum: f64 =
            gt.n.iter()
                .zip(&gt.proportion)
                .map(|(&n, &p)| n as f64 * p)
                .sum();
        assert!((sum - (1.0 - gt.p0)).abs() < 1e-9);
    }

    #[test]
    fn single_distinct_count_is_nan() {
        let gt = good_turing(&[5, 5, 5, 5, 5], 1.96);
        assert_eq!(gt.count, [5]);
        assert!(gt.proportion[0].is_nan());
        assert_eq!(gt.p0, 0.0);
    }

    #[test]
    fn zeros_feed_unseen_mass() {
        let gt = good_turing(&[0, 0, 0, 1, 1, 2, 3], 1.96);
        assert_eq!(gt.n0, 3);
        assert_eq!(gt.count, [1, 2, 3]);
        assert!((gt.p0 - 0.285714285714).abs() < 1e-9);
    }

    // Reference (r, N_r) tabulation independent of the max-vs-length branch,
    // so both the dense and sorted-unique paths can be checked against it.
    fn ref_tally(x: &[i64]) -> (i64, Vec<i64>, Vec<i64>) {
        use std::collections::BTreeMap;
        let mut m: BTreeMap<i64, i64> = BTreeMap::new();
        for &v in x {
            *m.entry(v).or_insert(0) += 1;
        }
        let n0 = m.remove(&0).unwrap_or(0);
        let (count, nr) = m.into_iter().unzip();
        (n0, count, nr)
    }

    // The dense tabulate (max_x < n) and the sorted-unique collector
    // (max_x >= n) must yield byte-identical (r, N_r) pairs; only the
    // collection strategy differs.
    #[test]
    fn both_tabulation_paths_agree() {
        // Dense path: max 3 < len 6.
        let dense = [1i64, 1, 2, 3, 3, 3];
        let gt = good_turing(&dense, 1.96);
        let (n0, count, nr) = ref_tally(&dense);
        assert_eq!(gt.n0, n0);
        assert_eq!(gt.count, count);
        assert_eq!(gt.n, nr);

        // Sorted path: max 2e9 >= len 5.
        let sparse = [2_000_000_000i64, 2_000_000_000, 1_000_000_000, 1, 2];
        let gt = good_turing(&sparse, 1.96);
        let (n0, count, nr) = ref_tally(&sparse);
        assert_eq!(gt.n0, n0);
        assert_eq!(gt.count, count);
        assert_eq!(gt.n, nr);
    }

    // A dense table sized to max(x) would allocate O(max_x) i64 — a single
    // ~2e9 count means ~16 GB and an OOM SIGKILL. The sorted-unique branch
    // makes this input finish in microseconds with the correct, finite SGT
    // result. Expected values are our i64-exact SGT output; edgeR reports
    // P0 = 1.4183739e-9 here because its C core sums the total sample size
    // (5_000_000_003) in a 32-bit int that wraps to 705_032_707 — a HELD
    // divergence in the > 2^31 total regime, where our i64 stays correct.
    #[test]
    fn large_count_column_no_oom() {
        let gt = good_turing(&[2_000_000_000, 2_000_000_000, 1_000_000_000, 1, 2], 1.96);
        assert_eq!(gt.n0, 0);
        assert_eq!(gt.count, [1, 2, 1_000_000_000, 2_000_000_000]);
        assert_eq!(gt.n, [1, 1, 1, 2]);
        assert_eq!(gt.p0, 1.0 / 5_000_000_003.0);
        let expect = [
            2.799_115_569_498_737e-10,
            4.869_217_280_128_239e-10,
            1.999_999_998_454_309_6e-1,
            3.999_999_995_938_678_4e-1,
        ];
        for (got, want) in gt.proportion.iter().zip(&expect) {
            assert!(
                (got - want).abs() <= want.abs() * 1e-12,
                "got {got} want {want}"
            );
        }
        let mass: f64 =
            gt.n.iter()
                .zip(&gt.proportion)
                .map(|(&n, &p)| n as f64 * p)
                .sum();
        assert!((mass - (1.0 - gt.p0)).abs() < 1e-15);
    }

    // A 1e10 count exceeds R's 32-bit integer domain entirely (as.integer(x)
    // returns NA there), so edgeR cannot process it at all — but the dense
    // table would try to allocate ~80 GB. The sorted-unique branch handles it
    // exactly; values are self-consistent SGT (proportions renormalize to
    // 1 - P0).
    #[test]
    fn count_beyond_r_int_domain_no_oom() {
        let gt = good_turing(&[10_000_000_000, 5_000_000_000, 3, 1, 1], 1.96);
        assert_eq!(gt.n0, 0);
        assert_eq!(gt.count, [1, 3, 5_000_000_000, 10_000_000_000]);
        assert_eq!(gt.n, [2, 1, 1, 1]);
        assert_eq!(gt.p0, 2.0 / 15_000_000_005.0);
        let expect = [
            9.104_490_002_888_423e-11,
            2.276_166_499_496_975_2e-10,
            3.333_333_331_623_107_3e-1,
            6.666_666_662_946_495e-1,
        ];
        for (got, want) in gt.proportion.iter().zip(&expect) {
            assert!(
                (got - want).abs() <= want.abs() * 1e-12,
                "got {got} want {want}"
            );
        }
        let mass: f64 =
            gt.n.iter()
                .zip(&gt.proportion)
                .map(|(&n, &p)| n as f64 * p)
                .sum();
        assert!((mass - (1.0 - gt.p0)).abs() < 1e-15);
    }
}
