# Efficient $H_0$ -based bitvector compression

This repository contains code related to our project to further research related to rrr_vectors from https://github.com/simongog/sdsl-lite

## Repo contents

### $H_0$ bitvector implementations
* `h0_bv.hpp`: implementation for decoding blocks using balansed recurrence.
* `h0_it.hpp`: implementation for decoding blocks iteratively in 8 block chunks with $b = 64$.
* `h0_63.hpp`: implementation for decoding blocks iteratively in 8 block chunks with $b = 63$.
* `h0_gap.hpp`: sparsified lookup table implementation.
* `weightedDeBrujin.py`: generates headers for our weighed deBrujin sequence based implementation.

### Hybrid bitvector implementations
* `hyb_vanilla.hpp`: same as the hybrid bitvector in the SDSL.
* `hyb_it.hpp`: hybrid bitvector implementation with option of encoding using the `h0_it` block implementation.
* `hyb_256.hpp`: **broken** hybrid bitvector implementation wth optino of encoding using the SDSL rrr_bitvector with $b = 256$.

### Benchmarking code
* `rrr_time_and_space.cpp`: Code to run benchmarks. Adapted from SDSL.
* `tester.sh`: Compiles and runs all test, generates test outputs in `res` folder.
* `Makefile`: Contains rules for creating all the benchmakring binaries.
* `results.ipynb`: Contains code for generating plots and exploring results.

Besides this there is debugging code, that is probably best left untouched untill it gets removed.

## Running

https://github.com/simongog/sdsl-lite should be installed on the system used for testing. The default install location is set up in the `Makefile`. 
If SDSL is not installed in the default location, the paths used in `ISDSL` and `SDSL_INCLUDE` should be updated.

With sdsl installed, tests can be run by calling `tester.sh <path/to/data/files>` in the repository root.

Results will be stored in the `res` folder.

Test data can be retrieved from:
* http://people.eng.unimelb.edu.au/sgog/data/WT-WEB-1GB.gz
* http://people.eng.unimelb.edu.au/sgog/data/WT-DNA-1GB.gz
* https://zenodo.org/records/11031752?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjZhN2Q4YTNmLTBhZmMtNDlmOS1iZjkwLWZlYmZjMjkzYzAzMSIsImRhdGEiOnt9LCJyYW5kb20iOiJlODNiYzkyYmQ4ZjRkZDYyMDcwZTFjZDRmODBmNWZmYiJ9.NaLDlFP16W_BV2jyMfenI42VU9SJLSEMELjaSxID7RfEMV1yM8Bhhmauf-6WvsEgmAlx6Bsx4BuRWG0rpTfiYA

