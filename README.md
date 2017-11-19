# metannot
Multithreaded wavelet trie construction library

Prerequisites
- C++ 11
- boost
- sdsl-lite

Usage
```
metannot NUM_BLOCKS CHECK_VECTOR INPUT_FILE1 [INPUT_FILE2 [...]] OUTPUT_FILE
```
Setting NUM_BLOCKS to 0 compresses the vector in a single block. CHECK_VECTOR is an indicator for whether the vector should be checked for successful compression after construction.

Input files are text files with one set of comma-separated numbers per line. Each line represents an input vector, while the numbers in a line indicate which bits in that line are set to 1.

If the `NJOBS=N` environment variable is set, then `N` threads will be used.
