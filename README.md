[![Crates.io](https://img.shields.io/crates/v/genomap)](https://crates.io/crates/genomap) [![Crates.io](https://img.shields.io/crates/d/genomap)](https://crates.io/crates/genomap) [![docs](https://docs.rs/genomap/badge.svg)](https://docs.rs/genomap) [![Rust CI](https://github.com/vsbuffalo/genomap/actions/workflows/rust.yml/badge.svg)](https://github.com/vsbuffalo/genomap/actions) [![DOI](https://zenodo.org/badge/754408774.svg)](https://zenodo.org/doi/10.5281/zenodo.10653718)



# A simple Rust library for storing data indexed by a chromosome name

genomap is a small library for storing a key-value map between chromosome
names and some generic data in a `GenomeMap`. Since in nearly every case we
want chromosomes to be sorted by their names, `GenomeMap` maintains an internal
sorted set of keys. `GenomeMap` uses a specialized chromosome name sorting
function that should properly sort autosomes, sex chromosomes, handle
Drosophila Chromosome names (e.g. 2L and 2R), etc. Please file [a GitHub
issue](http://github.com/vsbuffalo/genomemap/issues) if the sort order is not
as you'd anticipate.

Internally, the data stored in a `genomap::GenomeMap<T>` is in a `Vec<T>`, and
the type maintains a sorted list of chromosome names, and a forward and reverse
lookup table that associated the position in the `Vec` to the chromosome's name.

Below is a code example:


```rust
use genomap::GenomeMap;

let mut sm: GenomeMap<i32> = GenomeMap::new();
sm.insert("chr1", 1).unwrap();
sm.insert("chr2", 2).unwrap();

// get a reference to a value by name
println!("{:?}", sm.get("chr1"));

// iterate through name/values
for (name, value) in sm.iter() {
   println!("{} -> {}", name, value);
}

// get the index for a chromosome name
let index = sm.get_index_by_name("chr1").unwrap();
assert_eq!(index, 0);

// get a name by index
assert_eq!(sm.get_name_by_index(index).unwrap(), "chr1");
```

In Rust, working with non-`Copy`able types, such as a `String` chromosome name
key, can necessitate generic lifetime annotations. This can clutter code and
increase complexity significantly. To prevent this, `genomap` has *O(1)* access
by a `usize` index, so a chromosome name index can be stored in `Struct`s
rather than the `String` key.

## Performance

Multiple creation and access benchmarks are available in
`benches/comparison.rs`. Here is a small highlight of a sample of benchmarks.
For creation time, `GenomeMap` is about 20% slower. But this is incurred once
(and the absolute scale is insignificant).

| Data structure | Time      | Factor |
|----------------|-----------|--------|
| FnvHashMap     | 28.217 µs | 1.000  |
| IndexMap       | 29.844 µs | 1.058  |
| BTreeMap       | 29.401 µs | 1.041  |
| HashMap        | 29.420 µs | 1.043  |
| GenomeMap      | 33.913 µs | 1.202  |


`GenomeMap` has the second fastest sorted access times (it uses `FnvHashMap`'s
hasher internally, but there's one additional constant lookup time operation).

| Data structure | Time      | Factor |
|----------------|-----------|--------|
| FnvHashMap     | 68.555 ns | 1.00   |
| GenomeMap      | 198.55 ns | 2.89   |
| IndexMap       | 237.47 ns | 3.46   |
| HashMap        | 336.32 ns | 4.91   |
| BTreeMap       | 567.95 ns | 8.28   | 

