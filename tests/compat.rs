use std::collections::HashMap;
use std::process::Command;

fn ours() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_BIN_EXE_rsomics-edger-goodturing"))
}

fn golden(n: &str) -> String {
    format!("{}/tests/golden/{}", env!("CARGO_MANIFEST_DIR"), n)
}

fn run_ours(extra: &[&str]) -> String {
    let out = Command::new(ours())
        .arg(golden("counts.tsv"))
        .args(extra)
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "ours failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

fn parse_num(s: &str) -> f64 {
    if s == "NaN" {
        f64::NAN
    } else {
        s.parse().unwrap()
    }
}

/// Compare two TSVs cell-by-cell, keyed by (first column, header column), so a
/// column-reorder never causes a false diff. NaN matches NaN.
fn diff_matrix(mine: &str, theirs: &str, eps: f64) {
    let a = to_map(mine);
    let b = to_map(theirs);
    assert_eq!(a.len(), b.len(), "row count mismatch");
    let mut max_dev = 0.0f64;
    for (key, &x) in &a {
        let y = *b
            .get(key)
            .unwrap_or_else(|| panic!("cell {key:?} missing in oracle"));
        if x.is_nan() || y.is_nan() {
            assert!(
                x.is_nan() && y.is_nan(),
                "cell {key:?}: ours={x} oracle={y}"
            );
            continue;
        }
        let dev = (x - y).abs() / (1.0 + y.abs());
        max_dev = max_dev.max(dev);
        assert!(
            dev < eps,
            "cell {key:?}: ours={x} oracle={y} reldev={dev:e}"
        );
    }
    eprintln!("max relative deviation = {max_dev:e}");
}

fn to_map(s: &str) -> HashMap<(String, String), f64> {
    let mut lines = s.trim().lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').skip(1).collect();
    let mut out = HashMap::new();
    for line in lines {
        let mut f = line.split('\t');
        let row = f.next().unwrap().to_string();
        for (h, v) in header.iter().zip(f) {
            out.insert((row.clone(), h.to_string()), parse_num(v));
        }
    }
    out
}

/// The `gt` per-distinct-count table is keyed by (sample, r) so its float
/// columns compare independently of row order.
fn diff_gt(mine: &str, theirs: &str, eps: f64) {
    let a = gt_map(mine);
    let b = gt_map(theirs);
    assert_eq!(a.len(), b.len(), "gt row count mismatch");
    let mut max_dev = 0.0f64;
    for (key, va) in &a {
        let vb = b.get(key).unwrap_or_else(|| panic!("gt {key:?} missing"));
        for (col, (&x, &y)) in va.iter().zip(vb).enumerate() {
            if x.is_nan() || y.is_nan() {
                assert!(x.is_nan() && y.is_nan(), "gt {key:?} col {col}");
                continue;
            }
            let dev = (x - y).abs() / (1.0 + y.abs());
            max_dev = max_dev.max(dev);
            assert!(
                dev < eps,
                "gt {key:?} col {col}: ours={x} oracle={y} reldev={dev:e}"
            );
        }
    }
    eprintln!("gt max relative deviation = {max_dev:e}");
}

fn gt_map(s: &str) -> HashMap<(String, String), Vec<f64>> {
    let mut lines = s.trim().lines();
    lines.next();
    let mut out = HashMap::new();
    for line in lines {
        let f: Vec<&str> = line.split('\t').collect();
        let key = (f[0].to_string(), f[1].to_string());
        out.insert(
            key,
            vec![
                parse_num(f[2]),
                parse_num(f[3]),
                parse_num(f[4]),
                parse_num(f[5]),
            ],
        );
    }
    out
}

#[test]
fn proportions_match_golden() {
    diff_matrix(
        &run_ours(&["--proportions"]),
        &std::fs::read_to_string(golden("golden_prop.tsv")).unwrap(),
        1e-9,
    );
}

#[test]
fn per_count_matches_golden() {
    diff_gt(
        &run_ours(&[]),
        &std::fs::read_to_string(golden("golden_gt.tsv")).unwrap(),
        1e-9,
    );
}

// Live differential vs edgeR. Loud-skips when no edgeR Rscript is present.
#[test]
fn matches_edger_oracle_live() {
    let Some(rscript) = rscript() else {
        eprintln!("SKIP matches_edger_oracle_live: no edgeR Rscript");
        return;
    };
    let oracle = format!("{}/tests/goodturing_oracle.R", env!("CARGO_MANIFEST_DIR"));

    let prop = run_oracle(&rscript, &oracle, "prop", "1.96");
    diff_matrix(&run_ours(&["--proportions"]), &prop, 1e-9);

    let gt = run_oracle(&rscript, &oracle, "gt", "1.96");
    diff_gt(&run_ours(&[]), &gt, 1e-9);

    let gt3 = run_oracle(&rscript, &oracle, "gt", "3.0");
    diff_gt(&run_ours(&["--conf", "3.0"]), &gt3, 1e-9);
}

fn run_oracle(rscript: &[String], script: &str, mode: &str, conf: &str) -> String {
    let out = Command::new(&rscript[0])
        .args(&rscript[1..])
        .arg(script)
        .arg(golden("counts.tsv"))
        .arg(mode)
        .arg(conf)
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "oracle failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

fn rscript() -> Option<Vec<String>> {
    let home = std::env::var("HOME").unwrap_or_default();
    let candidates: Vec<Vec<String>> = vec![
        vec![format!("{home}/miniconda3/envs/r-bioc/bin/Rscript")],
        vec![
            "conda".into(),
            "run".into(),
            "-n".into(),
            "r-bioc".into(),
            "Rscript".into(),
        ],
        vec!["Rscript".into()],
    ];
    for c in candidates {
        let ok = Command::new(&c[0])
            .args(&c[1..])
            .arg("-e")
            .arg("suppressMessages(library(edgeR)); cat('ok')")
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false);
        if ok {
            return Some(c);
        }
    }
    None
}
