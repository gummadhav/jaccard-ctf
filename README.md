## Build instructions

This tool uses Cyclops Tensor Framework (CTF) routines. CTF is a distributed memory library that delivers routines for contraction and summation of sparse and dense tensors.

Download and build CTF: https://github.com/cyclops-community/ctf. 

Easiest way to do above is to install CTF with `./configure --no-dynamic && make` install globally with `sudo make install`, but local build would also work fine.

To compile this tool provide include and lib path to CTF in config.mk and run make.

## Run instructions

To test with default parameters, just run `./jaccard` or `mpirun -np 4 ./jaccard`.

To run, use e.g. `./jaccard -m 4000 -n 100 -p .01 -nbatch 10`, which would generate a 4000-by-100 k-mer bit matrix with 1% nonzeros, then compute a 100-by-100 similarity matrix by accumulating batches of 400 rows at a time.

To compute Jaccard index using k-mer files, run

`mpirun -np 4 ./jaccard -f "Path to directory that contains k-mer files" -lfile "List of k-mer files. Each file contains a data sample/read (One column of the input matrix, n)" -m "number of rows" -n "number of columns" -nbatch "number of batches"`
