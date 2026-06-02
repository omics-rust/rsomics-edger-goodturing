use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-edger-goodturing");
    let counts = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden/counts.tsv");
    c.bench_function("goodturing proportions golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .arg(counts.to_str().unwrap())
                .arg("--proportions")
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench);
criterion_main!(benches);
