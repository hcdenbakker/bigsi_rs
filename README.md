# bigsi_rs
Rust in-memory implementation of a BIGSI-like data structure

Rust in-memory implementation of a BIGSI-like data structure (see https://www.nature.com/articles/s41587-018-0010-1).
Comparable to a bloom filter; where a bloom filter tells if an element belongs to a single previously indexed set, a 
BIGSI-like data structure efficiently tells if an elements is a member of multiple query sets. Parameters (in particular
the index size and the number of hashes) should be chosen to assure a (down-stream) application permissable false positive probabilty. My strategy is to base this on the largest set to be indexed and use the formula of Goel and Gupta 2010 (https://web.stanford.edu/~ashishg/papers/inverted.pdf) to calculate the false postive probability for this set, and hence the maximum false positive probability for the index.     
  
# Example
```
use bigsi_rs;


//create a new index of size 250,000, 10 accessions and 3 hash functions
let mut new_filter = bigsi_rs::Bigsi::new(250000, 10, 3);

//insert words in index for accessions 0, 3 and 7        
new_filter.insert(0, "ATGT");
new_filter.insert(3, "ATGT");
new_filter.insert(7, "ATGT");

//shrink uninformative elements in index
new_filter.slim();

assert_eq!(new_filter.get("ATGT").len(), 3 as usize);
assert_eq!(new_filter.get("ATGC").len(), 0 as usize);
```
