use criterion::{black_box, criterion_group, criterion_main, Criterion};
use fnv::{FnvBuildHasher, FnvHashMap};
use indexmap::IndexMap;
use rand::{thread_rng, Rng};
use std::collections::{BTreeMap, HashMap};

use genomap::GenomeMap;

pub fn chrom(num: usize) -> String {
    format!("chr{}", num)
}

pub fn random_chrom(nchrom: usize) -> String {
    let mut rng = thread_rng();
    chrom(rng.gen_range(1..=nchrom))
}

pub fn random_vec(nsize: usize) -> Vec<f64> {
    let mut rng = thread_rng();
    (0..nsize).map(|_| rng.gen::<f64>()).collect()
}

pub fn criterion_benchmark_creation(c: &mut Criterion) {
    let mut creation = c.benchmark_group("Creation");

    // GenomeMap creation benchmark
    creation.bench_function("GenomeMap creation 22", |b| {
        b.iter(|| {
            let mut gm = GenomeMap::new();
            for i in 1..=black_box(22) {
                gm.insert(&chrom(i), random_vec(black_box(100))).unwrap();
            }
        })
    });

    // HashMap creation benchmark
    creation.bench_function("HashMap creation 22", |b| {
        b.iter(|| {
            let mut hm = HashMap::new();
            for i in 1..=black_box(22) {
                hm.insert(chrom(i), random_vec(black_box(100)));
            }
        })
    });

    // FnvHashMap creation benchmark
    creation.bench_function("FnvHashMap creation 22", |b| {
        b.iter(|| {
            let mut hm: FnvHashMap<String, Vec<f64>> = FnvHashMap::default();
            for i in 1..=black_box(22) {
                hm.insert(chrom(i), random_vec(black_box(100)));
            }
        })
    });

    // IndexMap creation benchmark
    creation.bench_function("IndexMap creation 22", |b| {
        b.iter(|| {
            let mut hm: IndexMap<String, Vec<f64>> = IndexMap::new();
            for i in 1..=black_box(22) {
                hm.insert(chrom(i), random_vec(black_box(100)));
            }
        })
    });

    // BTreeMap creation benchmark
    creation.bench_function("BTreeMap creation 22", |b| {
        b.iter(|| {
            let mut hm: BTreeMap<String, Vec<f64>> = BTreeMap::new();
            for i in 1..=black_box(22) {
                hm.insert(chrom(i), random_vec(black_box(100)));
            }
        })
    });
}

pub fn criterion_benchmark_access(c: &mut Criterion) {
    let mut access = c.benchmark_group("Access");

    // GenomeMap access benchmark
    access.bench_function("GenomeMap access 22", |b| {
        let mut gm = GenomeMap::new();
        for i in 1..=22 {
            gm.insert(&chrom(i), random_vec(100)).unwrap();
        }
        b.iter(|| {
            gm.get(&random_chrom(22)).unwrap();
        })
    });

    // HashMap access benchmark
    access.bench_function("HashMap access 22", |b| {
        let mut hm = HashMap::new();
        for i in 1..=22 {
            hm.insert(chrom(i), random_vec(100));
        }
        b.iter(|| {
            hm.get(&random_chrom(22)).unwrap();
        })
    });

    // FnvHashMap access benchmark
    access.bench_function("FnvHashMap access 22", |b| {
        let mut hm: FnvHashMap<String, Vec<f64>> = FnvHashMap::default();
        for i in 1..=22 {
            hm.insert(chrom(i), random_vec(100));
        }
        b.iter(|| {
            hm.get(&random_chrom(22)).unwrap();
        })
    });

    // IndexMap access benchmark
    access.bench_function("IndexMap access 22", |b| {
        let mut hm: IndexMap<String, Vec<f64>> = IndexMap::new();
        for i in 1..=22 {
            hm.insert(chrom(i), random_vec(100));
        }
        b.iter(|| {
            hm.get(&random_chrom(22)).unwrap();
        })
    });

    // BTreeMap access benchmark
    access.bench_function("BTreeMap access 22", |b| {
        let mut hm: BTreeMap<String, Vec<f64>> = BTreeMap::new();
        for i in 1..=22 {
            hm.insert(chrom(i), random_vec(100));
        }
        b.iter(|| {
            hm.get(&random_chrom(22)).unwrap();
        })
    });
}

pub fn criterion_benchmark_sorted(c: &mut Criterion) {
    let mut sorted = c.benchmark_group("Sorted");

    // GenomeMap access benchmark with sorted keys
    sorted.bench_function("GenomeMap sorted access 22", |b| {
        let mut gm = GenomeMap::new();
        for i in 1..=22 {
            gm.insert(&chrom(i), random_vec(100)).unwrap();
        }
        let sorted_keys: Vec<_> = gm.names();
        b.iter(|| {
            for key in sorted_keys.iter() {
                black_box(gm.get(key).unwrap());
            }
        })
    });

    // HashMap access benchmark with sorted keys
    sorted.bench_function("HashMap sorted access 22", |b| {
        let mut hm = HashMap::new();
        for i in 1..=22 {
            hm.insert(chrom(i), random_vec(100));
        }
        let sorted_keys: Vec<_> = hm.keys().cloned().collect();
        b.iter(|| {
            for key in sorted_keys.iter() {
                black_box(hm.get(key).unwrap());
            }
        })
    });

    // FnvHashMap access benchmark with sorted keys
    sorted.bench_function("FnvHashMap sorted access 22", |b| {
        let mut hm: FnvHashMap<String, Vec<f64>> = FnvHashMap::default();
        for i in 1..=22 {
            hm.insert(chrom(i), random_vec(100));
        }
        let sorted_keys: Vec<_> = hm.keys().cloned().collect();
        b.iter(|| {
            for key in sorted_keys.iter() {
                black_box(hm.get(key).unwrap());
            }
        })
    });

    // IndexMap access benchmark with sorted keys
    sorted.bench_function("IndexMap sorted access 22", |b| {
        let mut hm: IndexMap<String, Vec<f64>> = IndexMap::new();
        for i in 1..=22 {
            hm.insert(chrom(i), random_vec(100));
        }
        // IndexMap maintains insertion order, so if insertion was sorted, this step might be unnecessary
        // But to ensure it, we explicitly sort the keys here as well
        let sorted_keys: Vec<_> = hm.keys().cloned().collect();
        b.iter(|| {
            for key in sorted_keys.iter() {
                black_box(hm.get(key).unwrap());
            }
        })
    });

    // BTreeMap access benchmark with sorted keys
    sorted.bench_function("BTreeMap sorted access 22", |b| {
        let mut hm: BTreeMap<String, Vec<f64>> = BTreeMap::new();
        for i in 1..=22 {
            hm.insert(chrom(i), random_vec(100));
        }

        let keys: Vec<_> = hm.keys().cloned().collect();
        b.iter(|| {
            for key in &keys {
                black_box(hm.get(key).unwrap());
            }
        })
    });
}

criterion_group!(
    benches,
    criterion_benchmark_creation,
    criterion_benchmark_access,
    criterion_benchmark_sorted
);
criterion_main!(benches);
