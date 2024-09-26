//! A fast, ordered map-like data structure for storing genomic data associated with chromosomes.
//!
//! The only data structure is a [`GenomeMap`] that is like an [`HashMap`], but with
//! the chromosome name keys automatically ordered according to the usual biological order
//! (see below for details). Unlike `HashMap`, `GenomeMap` will also not let you clobber
//! existing entries with [`GenomeMap.insert()`] by raising an error. Here is a basic example:
//!
//! ```
//! use genomap::GenomeMap;
//!
//! let mut sm: GenomeMap<i32> = GenomeMap::new();
//! sm.insert("chr1", 1).unwrap();
//! sm.insert("chr2", 2).unwrap();
//! println!("{:?}", sm.get("chr1"));
//! ```
//!
//! In Rust, since working with non-`Copy`able types (such as a `String` chromosome name)
//! can necessitate generic lifetime annotations or cloning. This can clutter code
//! and increase complexity, or decrease performance. To alleviate this, `genomap` supports *O(1)* access
//! by a `usize` index corresponding to the key at that index in the sorted names. This
//! index can be stored in a `Struct`, which prevents (1) cloning lots of `String` names and (2)
//! managing the lifetimes to `&str` in types.
//!
//! ```
//! use genomap::GenomeMap;
//!
//! let mut sm: GenomeMap<i32> = GenomeMap::new();
//! sm.insert("chr1", 1).unwrap();
//! sm.insert("chr2", 2).unwrap();
//!
//! for (name, value) in sm.iter() {
//!    println!("{} -> {}", name, value);
//! }
//!
//! let index = sm.get_index_by_name("chr1").unwrap();
//! assert_eq!(index, 0);
//!
//! assert_eq!(sm.get_name_by_index(index).unwrap(), "chr1");
//!
//! ```
//!
//! ## Chromosome Name Ordering
//!
//! This sorts the names (e.g. chromosome or contig name) by doing the following:
//! 1. Remove 'chr' if present.
//! 2. See if the remaining bit is a number - if so, that comes first. To handle *Drosophila*
//!    chromosomes, it has a special rule to sort 2L and 2R into the integer autosomes.
//! 3. Then, letters. But: order common sex chromosome names (X, Y, M, Z, W, O), then mitochondria,
//!    then everything else.
//!
//! If the ordering created by this system is not what you'd expect for your organism, please
//! file an issue on GitHub: <http://github.com/vsbuffalo/genomemap/issues>
//!
//! [`GenomeMap.insert()`]: crate::GenomeMap::insert

use fnv::FnvBuildHasher;
use std::{cmp::Ordering, collections::HashMap};
use thiserror::Error;

/// [`GenomeMap`] is an ordered, map-like data structure for storing genomic data associated with
/// chromosomes.
///
/// [`GenomeMap`] maintains an key-value mapping, much like an
/// `[IndexMap](https://docs.rs/indexmap/1.9.3/indexmap/)`. However, while an `IndexMap` is ordered
/// by insertion, the data we store in genomes should be ordered according to chromosome or contig
/// name order.
///
#[derive(Clone, Debug)]
pub struct GenomeMap<T> {
    // relation between name and insertion order index
    lookup: HashMap<String, usize, FnvBuildHasher>,
    // relation between insertion order index and name
    reverse_lookup: HashMap<usize, String, FnvBuildHasher>,
    // values
    values: Vec<T>,
    // the always-ordered set of keys; their position is
    // the external indices.
    pub sorted_keys: Vec<String>,
}

// Note: cannot use derive macro, as this requires T implements default
impl<T> Default for GenomeMap<T> {
    fn default() -> Self {
        Self {
            lookup: HashMap::default(),
            reverse_lookup: HashMap::default(),
            values: Vec::default(),
            sorted_keys: Vec::with_capacity(50),
        }
    }
}

#[derive(Debug, Error)]
pub enum GenomeMapError {
    #[error("GenomeMap.insert() error: already contains the sequence '{0}'")]
    SequenceMapContainsSeqname(String),
}

impl<W, T> PartialEq<GenomeMap<W>> for GenomeMap<T>
where
    T: PartialEq<W>,
    W: IntoIterator,
    T: IntoIterator,
{
    fn eq(&self, other: &GenomeMap<W>) -> bool {
        self.iter()
            .zip(other.iter())
            .all(|((k1, v1), (k2, v2))| k1 == k2 && v1 == v2)
    }
}

impl<T> GenomeMap<T> {
    /// Create a new [`GenomeMap`].
    pub fn new() -> Self {
        GenomeMap::default()
    }

    /// Retrieve a reference to the value stored for `name`, if it is present, else `None`.
    pub fn get(&self, name: &str) -> Option<&T> {
        self.lookup.get(name).map(|idx| &self.values[*idx])
    }

    /// Retrieve a *mutable* reference to the value stored for `name`, if it is present, else `None`.
    pub fn get_mut(&mut self, name: &str) -> Option<&mut T> {
        self.lookup.get_mut(name).map(|idx| &mut self.values[*idx])
    }

    /// Retrieve a reference to the value stored for the specified `index`.
    pub fn get_by_index(&self, index: usize) -> Option<&T> {
        // get the name from the *external* index, i.e. the order of the names
        let name = self.get_name_by_index(index);
        name.and_then(|key| self.get(key))
    }

    /// Return a reference (i.e. `&str`) to the name (e.g. chromosome or contig name)
    /// for the specified index. If not present, will return `None`. This is *O(1)*.
    pub fn get_name_by_index(&self, index: usize) -> Option<&str> {
        self.sorted_keys.get(index).as_ref().map(|x| x.as_str())
    }

    /// Return the index for a particular contig, by doing binary search.
    /// If the name is not found, returns `None`. This is *O(log(n))*
    /// where *n* is the number of elements.
    pub fn get_index_by_name(&self, name: &str) -> Option<usize> {
        match self
            .sorted_keys
            .binary_search_by(|probe| chromosome_probe(probe, name))
        {
            Ok(index) => Some(index),
            Err(_) => None,
        }
    }

    /// Get the ordered names. Their order corresponds to the indices.
    pub fn names(&self) -> Vec<String> {
        self.sorted_keys.clone()
    }

    /// Get a slice of references to the ordered names.
    pub fn names_ref(&self) -> &[String] {
        &self.sorted_keys
    }

    /// Get the length of the [`GenomeMap`].
    pub fn len(&self) -> usize {
        let len = self.lookup.len();
        assert_eq!(len, self.values.len());
        len
    }

    /// Is the container empty?
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Get the indices of the names (this is just `0..GenomeMap.len()`).
    pub fn indices(&self) -> std::ops::Range<usize> {
        std::ops::Range {
            start: 0,
            end: self.len(),
        }
    }

    /// Remove an entry from the GenomeMap by its key (chromosome name).
    /// Returns the value if the key was present in the map, or None if it wasn't.
    pub fn remove(&mut self, name: &str) -> Option<T> {
        if let Some(&index) = self.lookup.get(name) {
            // Remove from lookup
            self.lookup.remove(name);

            // Remove from sorted_keys
            if let Ok(sorted_index) = self.sorted_keys.binary_search(&name.to_string()) {
                self.sorted_keys.remove(sorted_index);
            }

            // Remove the value
            let value = self.values.remove(index);

            // Update indices for all elements after the removed one
            for i in index..self.values.len() {
                if let Some(chrom_name) = self.reverse_lookup.get(&(i + 1)) {
                    self.lookup.insert(chrom_name.clone(), i);
                    self.reverse_lookup.insert(i, chrom_name.clone());
                }
            }

            // Remove the last reverse_lookup entry
            self.reverse_lookup.remove(&self.values.len());

            Some(value)
        } else {
            None
        }
    }

    /// O(1) check if a particular name is in the container.
    pub fn contains(&self, name: &str) -> bool {
        self.lookup.contains_key(name)
    }

    /// Insert a value into the [`GenomeMap`] for the given name. Will error if an entry
    /// with this name already exists.
    pub fn insert(&mut self, name: &str, value: T) -> Result<usize, GenomeMapError> {
        if self.contains(name) {
            return Err(GenomeMapError::SequenceMapContainsSeqname(name.to_string()));
        }
        self.values.push(value);
        // this is the *internal* index, based on the position in the Vec
        let insertion_index = self.values.len() - 1;
        self.lookup.insert(name.to_string(), insertion_index);
        self.reverse_lookup
            .insert(insertion_index, name.to_string());

        // insert into sorted_keys and maintain the order: this is the *external* index
        let index = match self
            .sorted_keys
            .binary_search_by(|probe| chromosome_probe(probe, name))
        {
            Ok(_) => {
                panic!("An internal error was encountered. Please report an issue: https://github.com/vsbuffalo/genomemap/issues");
            }
            Err(pos) => {
                self.sorted_keys.insert(pos, name.to_string());
                pos
            }
        };
        Ok(index)
    }

    /// Return a mutable reference to a value if the name is found,
    /// otherwise return a new default value.
    pub fn entry_or_default(&mut self, name: &str) -> &mut T
    where
        T: Default,
    {
        let insertion_index = match self.lookup.get(name) {
            Some(&index) => index,
            None => {
                let new_entry = T::default();
                let _ = self.insert(name, new_entry).unwrap();
                // get the insertion_index in the lookup table
                *self.lookup.get(name).unwrap()
            }
        };

        &mut self.values[insertion_index]
    }

    /// Iterate over `(name, &value)` tuples, ordered by name.
    pub fn iter(&self) -> impl Iterator<Item = (&String, &T)> {
        self.sorted_keys
            .iter()
            .map(|name| (name, self.get(name).unwrap()))
    }

    /// Iterate over mutable values.
    pub fn values_mut(&mut self) -> std::slice::IterMut<'_, T> {
        self.values.iter_mut()
    }

    /// Get the values, ordered according to the sorted names order.
    pub fn values(&self) -> impl Iterator<Item = &T> {
        self.sorted_keys.iter().map(|name| self.get(name).unwrap())
    }

    /// Shrink all internal data structures to fit.
    pub fn shrink_to_fit(&mut self) {
        self.lookup.shrink_to_fit();
        self.reverse_lookup.shrink_to_fit();
        self.sorted_keys.shrink_to_fit();
    }
}

impl<T> IntoIterator for GenomeMap<T> {
    type Item = (String, T);
    type IntoIter = std::vec::IntoIter<Self::Item>;

    /// Make a consuming iterator over the `(String, T)` keys and values.
    fn into_iter(self) -> Self::IntoIter {
        let pairs = self
            .sorted_keys
            .into_iter()
            .zip(self.values)
            .collect::<Vec<_>>();
        pairs.into_iter()
    }
}

/// For `collect()`ing name-value tuples into a new [`GenomeMap`].
///
/// # Panics
/// This will panic if iterator processor encounters a sequence already processed
/// in the iterator. This is a protection against accidental clobbering. Cases when users would not
/// want this likely reflect code smell.
impl<'a, T> FromIterator<(String, T)> for GenomeMap<T>
where
    T: 'a,
{
    fn from_iter<I: IntoIterator<Item = (String, T)>>(iter: I) -> Self {
        let mut genomap = GenomeMap::new();

        for (name, item) in iter {
            let msg = format!(
                "Invalid iterator: contains multiple entries for '{}', which would be clobbered.",
                name
            );
            genomap.insert(&name, item).expect(&msg);
        }
        genomap
    }
}

/// A probe for `binary_search_by` that orders the names, according to typical chromosome or contig
/// naming order.
///
/// This sorts by doing the following:
/// 1. Remove 'chr' if there.
/// 2. See if the remaining bit is a number - if so, that comes first.
/// 3. Then, letters. But: order common sex chromosome names (X, Y, M, Z, W, O), then mitochondria,
///    then everything else.
pub fn chromosome_probe(name: &str, target: &str) -> Ordering {
    let extract_parts = |chr: &str| -> (u8, Option<u32>, String, u8) {
        let chr_trimmed = chr.trim_start_matches("chr");
        match chr_trimmed {
            "X" => (2, None, chr_trimmed.to_string(), 0),
            "Y" => (3, None, chr_trimmed.to_string(), 0),
            "M" | "MT" | "Mt" => (4, None, chr_trimmed.to_string(), 0),
            "Z" => (5, None, chr_trimmed.to_string(), 0),
            "W" => (6, None, chr_trimmed.to_string(), 0),
            "O" => (7, None, chr_trimmed.to_string(), 0),
            _ if chr_trimmed.ends_with('L') || chr_trimmed.ends_with('R') => {
                let (num, arm) = chr_trimmed.split_at(chr_trimmed.len() - 1);
                if let Ok(num) = num.parse::<u32>() {
                    let arm_priority = if arm == "L" { 1 } else { 2 }; // L before R
                    (1, Some(num), String::new(), arm_priority)
                } else {
                    (8, None, chr_trimmed.to_string(), 0)
                }
            }
            _ if chr_trimmed.parse::<u32>().is_ok() => {
                (1, Some(chr_trimmed.parse().unwrap()), String::new(), 0)
            }
            _ => (8, None, chr_trimmed.to_string(), 0),
        }
    };

    let name_parts = extract_parts(name);
    let target_parts = extract_parts(target);

    name_parts.cmp(&target_parts)
}

#[cfg(test)]
mod test {
    use super::GenomeMap;

    #[test]
    fn test_genomemap_new() {
        let mut sm: GenomeMap<i32> = GenomeMap::new();
        assert!(sm.is_empty());
        sm.insert("chr1", 1).unwrap();
        sm.insert("chr2", 2).unwrap();

        assert_eq!(sm.indices(), 0..2);

        assert_eq!(sm.len(), 2);
        assert!(!sm.is_empty());

        // check that it insertion order-invariant.
        let mut sm: GenomeMap<i32> = GenomeMap::new();
        sm.insert("chr2", 1).unwrap();
        sm.insert("chr1", 2).unwrap();

        assert_eq!(sm.indices(), 0..2);
    }

    #[test]
    fn test_genomemap_contains() {
        let mut sm: GenomeMap<i32> = GenomeMap::new();
        assert!(sm.is_empty());
        sm.insert("chr1", 1).unwrap();
        sm.insert("chr2", 2).unwrap();
        assert!(sm.contains("chr1"));
        assert!(sm.contains("chr2"));
        assert!(!sm.contains("chrX"));
    }

    #[test]
    fn test_genomemap_shrink() {
        let mut sm: GenomeMap<i32> = GenomeMap::new();
        sm.insert("chr1", 1).unwrap();
        sm.insert("chr2", 2).unwrap();
        sm.shrink_to_fit();
    }

    #[test]
    #[should_panic]
    fn test_genomemap_clobber() {
        let mut sm: GenomeMap<i32> = GenomeMap::new();
        sm.insert("chr1", 1).unwrap();
        sm.insert("chr1", 2).unwrap();
    }

    #[test]
    fn test_entry_or_default() {
        let mut sm: GenomeMap<Vec<i32>> = GenomeMap::new();
        let vec = sm.entry_or_default("chr2");
        assert!(vec.is_empty());
        vec.push(1);
        vec.push(2);
        let vec = sm.entry_or_default("chr1");
        vec.push(1);

        let vec_chr1 = sm.get("chr2");
        assert_eq!(*vec_chr1.unwrap(), vec![1, 2]);

        assert_eq!(sm.len(), 2);
    }

    #[test]
    fn test_genomemap_get_mut() {
        let mut sm: GenomeMap<i32> = GenomeMap::new();
        sm.insert("chr1", 1).unwrap();
        sm.insert("chr2", 2).unwrap();

        assert_eq!(*sm.get("chr1").unwrap(), 1);
        let val = sm.get_mut("chr1").unwrap();
        *val = 10;
        assert_eq!(*sm.get("chr1").unwrap(), 10);
    }

    #[test]
    fn test_genomemap_get() {
        let mut sm: GenomeMap<i32> = GenomeMap::new();
        sm.insert("chr1", 1).unwrap();
        sm.insert("chr2", 2).unwrap();

        assert_eq!(*sm.get("chr1").unwrap(), 1);
        assert_eq!(*sm.get("chr2").unwrap(), 2);
        assert_eq!(sm.get("chr4"), None);

        // again, check if insertion order-invariant
        let mut sm: GenomeMap<i32> = GenomeMap::new();
        sm.insert("chr2", 2).unwrap();
        sm.insert("chr1", 1).unwrap();

        assert_eq!(*sm.get("chr1").unwrap(), 1);
        assert_eq!(*sm.get("chr2").unwrap(), 2);
        assert_eq!(sm.get("chr4"), None);
    }

    #[test]
    fn test_genomemap_get_by_index() {
        let mut sm: GenomeMap<i32> = GenomeMap::new();
        sm.insert("chr1", 1).unwrap();
        sm.insert("chr2", 2).unwrap();

        let chr1_idx = sm.get_index_by_name("chr1");
        assert_eq!(chr1_idx.unwrap(), 0);

        let chr2_idx = sm.get_index_by_name("chr2");
        assert_eq!(chr2_idx.unwrap(), 1);

        let chr3_idx = sm.get_index_by_name("chr3");
        assert_eq!(chr3_idx, None);

        assert_eq!(*sm.get("chr1").unwrap(), 1);
        assert_eq!(*sm.get("chr2").unwrap(), 2);
        assert_eq!(sm.get("chr3"), None);

        assert_eq!(*sm.get_by_index(chr1_idx.unwrap()).unwrap(), 1);
        assert_eq!(*sm.get_by_index(chr2_idx.unwrap()).unwrap(), 2);
    }

    #[test]
    fn test_map_get_name_by_index_and_name() {
        let mut sm: GenomeMap<i32> = GenomeMap::new();
        sm.insert("chr1", 1).unwrap();
        sm.insert("chr2", 2).unwrap();

        let index = sm.get_index_by_name("chr1").unwrap();
        assert_eq!(*sm.get_by_index(index).unwrap(), 1);

        assert_eq!(sm.get_name_by_index(index).unwrap(), "chr1");
    }

    #[test]
    fn test_map_get_name_by_index() {
        let mut sm: GenomeMap<i32> = GenomeMap::new();
        sm.insert("chr1", 1).unwrap();
        sm.insert("chr2", 2).unwrap();

        assert_eq!(sm.get_name_by_index(0).unwrap(), "chr1");
        assert_eq!(sm.get_name_by_index(1).unwrap(), "chr2");

        // should be insertion order permutation invariant
        let mut sm: GenomeMap<i32> = GenomeMap::new();
        sm.insert("chr2", 2).unwrap();
        sm.insert("chr1", 1).unwrap();

        assert_eq!(sm.get_name_by_index(0).unwrap(), "chr1");
        assert_eq!(sm.get_name_by_index(1).unwrap(), "chr2");
    }

    #[test]
    fn test_sorting_human() {
        let mut sm: GenomeMap<i32> = GenomeMap::new();

        sm.insert("chr2", 2).unwrap();
        sm.insert("chr1", 1).unwrap();
        sm.insert("chr22", 3).unwrap();
        sm.insert("chrX", 4).unwrap();
        sm.insert("chrY", 5).unwrap();
        sm.insert("chrMt", 6).unwrap();

        let expected_order = vec!["chr1", "chr2", "chr22", "chrX", "chrY", "chrMt"];
        assert_eq!(sm.names(), expected_order);

        // for ensembl-style naming
        let mut sm: GenomeMap<i32> = GenomeMap::new();

        sm.insert("2", 2).unwrap();
        sm.insert("1", 1).unwrap();
        sm.insert("22", 3).unwrap();
        sm.insert("X", 4).unwrap();
        sm.insert("Y", 5).unwrap();
        sm.insert("Mt", 6).unwrap();

        let expected_order = vec!["1", "2", "22", "X", "Y", "Mt"];
        assert_eq!(sm.names(), expected_order);
    }

    #[test]
    fn test_sorting_drosophila() {
        let mut sm: GenomeMap<i32> = GenomeMap::new();

        sm.insert("chr2L", 2).unwrap();
        sm.insert("chr1", 1).unwrap();
        sm.insert("chr2R", 3).unwrap();
        sm.insert("chrX", 5).unwrap();
        sm.insert("chrY", 6).unwrap();
        sm.insert("chr4", 4).unwrap();

        let expected_order = vec!["chr1", "chr2L", "chr2R", "chr4", "chrX", "chrY"];
        assert_eq!(sm.names(), expected_order);
    }

    #[test]
    fn test_into_iter() {
        let mut sm: GenomeMap<i32> = GenomeMap::new();
        sm.insert("chr1", 1).unwrap();
        sm.insert("chr2", 2).unwrap();
        let mut iter = sm.into_iter();
        assert_eq!(iter.next(), Some(("chr1".to_string(), 1)));
        assert_eq!(iter.next(), Some(("chr2".to_string(), 2)));
        assert!(iter.next().is_none());
    }

    #[test]
    fn test_genomemap_remove() {
        let mut sm: GenomeMap<i32> = GenomeMap::new();
        sm.insert("chr1", 1).unwrap();
        sm.insert("chr2", 2).unwrap();
        sm.insert("chr3", 3).unwrap();
        sm.insert("chrX", 4).unwrap();
        sm.insert("chrY", 5).unwrap();

        // Test removal of existing entry
        assert_eq!(sm.remove("chr2"), Some(2));
        assert_eq!(sm.len(), 4);
        assert!(!sm.contains("chr2"));

        // Test removal of non-existing entry
        assert_eq!(sm.remove("chr4"), None);
        assert_eq!(sm.len(), 4);

        // Check if remaining entries are still valid and in correct order
        assert_eq!(sm.names(), vec!["chr1", "chr3", "chrX", "chrY"]);
        assert_eq!(*sm.get("chr1").unwrap(), 1);
        assert_eq!(*sm.get("chr3").unwrap(), 3);
        assert_eq!(*sm.get("chrX").unwrap(), 4);
        assert_eq!(*sm.get("chrY").unwrap(), 5);

        // Test removal of first entry
        assert_eq!(sm.remove("chr1"), Some(1));
        assert_eq!(sm.len(), 3);
        assert_eq!(sm.names(), vec!["chr3", "chrX", "chrY"]);

        // Test removal of last entry
        assert_eq!(sm.remove("chrY"), Some(5));
        assert_eq!(sm.len(), 2);
        assert_eq!(sm.names(), vec!["chr3", "chrX"]);

        // Check if indices are still correct
        assert_eq!(sm.get_index_by_name("chr3"), Some(0));
        assert_eq!(sm.get_index_by_name("chrX"), Some(1));
    }

    #[test]
    fn test_genomemap_remove_and_insert() {
        let mut sm: GenomeMap<i32> = GenomeMap::new();
        sm.insert("chr1", 1).unwrap();
        sm.insert("chr2", 2).unwrap();
        sm.insert("chr3", 3).unwrap();

        // Remove middle entry
        sm.remove("chr2");

        // Insert new entry
        sm.insert("chr4", 4).unwrap();

        // Check if order and values are correct
        assert_eq!(sm.names(), vec!["chr1", "chr3", "chr4"]);
        assert_eq!(*sm.get("chr1").unwrap(), 1);
        assert_eq!(*sm.get("chr3").unwrap(), 3);
        assert_eq!(*sm.get("chr4").unwrap(), 4);

        // Check indices
        assert_eq!(sm.get_index_by_name("chr1"), Some(0));
        assert_eq!(sm.get_index_by_name("chr3"), Some(1));
        assert_eq!(sm.get_index_by_name("chr4"), Some(2));
    }
}
