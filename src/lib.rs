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

    let mut tab = vec![0i64; (max_x + 1) as usize];
    for &v in x {
        tab[v as usize] += 1;
    }
    let n0 = tab[0];

    let mut count = Vec::new();
    let mut nr = Vec::new();
    for (r, &c) in tab.iter().enumerate().skip(1) {
        if c > 0 {
            count.push(r as i64);
            nr.push(c);
        }
    }

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
}
