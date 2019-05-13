#main compiler to use
CXX=mpicxx
#compiler flags (-g -O0 for debug, -O3 for optimization), generally need -fopenmp and -std=c++0x
CXXFLAGS=-std=c++0x -g -O0 -fopenmp
#path to CTF include directory prefixed with -I
INCLUDE_PATH=-I/home/edgar/work/ctf/include
#path to MPI/CTF/scalapack/lapack/blas library directories prefixed with -L
LIB_PATH = -L/home/edgar/work/ctf/scalapack/build/lib -L/opt/intel/mkl/lib/intel64/ -L/home/edgar/work/ctf/hptt/lib
#libraries to link (MPI/CTF/scalapack/lapack/blas) to prefixed with -l
LIB_FILES = -lctf -Wl,-Bstatic -lscalapack -Wl,-Bdynamic -Wl,-Bstatic -Wl,--start-group -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -Wl,-Bdynamic -Wl,-Bstatic -lhptt -Wl,-Bdynamic
LINKFLAGS = -lpthread -lm -ldl -lgfortran
LIBS = $(LIB_FILES) $(LINKFLAGS)
