## Build instructions

Build CTF, provide include and lib path to config.mk and run make.

Easiest way to do above is to install CTF with `./configure --no-dynamic && make` install globally with `sudo make install`, but local build would also work fine.

## Run instructions

To test with default parameters, just run `./jaccard` or `mpirun -np 4 ./jaccard`.

To run, use e.g. `./jaccard -m 4000 -n 100 -p .01 -nbatch 10`, which would generate a 4000-by-100 k-mer bit matrix with 1% nonzeros, then compute a 100-by-100 similarity matrix by accumulating batches of 400 rows at a time.

