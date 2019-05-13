#include <ctf.hpp>

using namespace CTF;

/**
 * \brief get number of set bits in uint-like type bitmask
 *        TODO: when possible, can replace this with hardware instructions
 * \param[in] x bitmask cotnaining some 0 and 1 bits
 * \return number of 1 bits in x
 */
template <typename bitmask>
uint32_t popcount(bitmask x){
  uint32_t p = 0;
  for (int i=0; i<sizeof(bitmask)*8; i++){
    p += (x>>i) & 1;
  }
  return p;
}

/**
 * \brief Generate an m-by-n bool sparse matrix, then, given bitmasks of len_bm bits,
 *        generate a random (m/len_bm)-by-n sparse CTF matrix,
 *        where every element is a bitmask for a subcolumn of length k that contains at least 1 bit
 * \param[in] m number of rows in overall bool matrix
 * \param[in] n number of columns in matrix
 * \param[in] p probability of any element being nonzero (even if =1, matrix may have some zeros)
 * \param[in] dw CTF World (MPI comm) on which the matrix should live
 * \return ceil(m/len_bm)-by-n sparse CTF matrix where a_ij is a bitmask encoding a subcolumn of len_bm bits
 */
template <typename bitmask>
Matrix<bitmask> generate_random_A(int64_t m, int64_t n, double p, World & dw){
  int len_bm = sizeof(bitmask)*8;
  int64_t mm = (m + len_bm - 1)/len_bm;
  // define CTF sparse matrix
  Matrix<bitmask> A(mm,n,SP,dw);
  // define vector of (int64_t,bitmask) pairs, storing entries of A indexed as row+col*nrow
  // TODO: filling a large vector can require overhead due to reallocation/copy, should predefine size
  std::vector<Pair<bitmask>> pairs;
  // define communicator containing each processor by itself
  World selfw(MPI_COMM_SELF);
  for (int64_t i=dw.rank; i<n; i+=dw.np){
    // define a bool CTF sparse vector locally, used only to generate a random sparse vector
    Vector<bool> v(m,SP,selfw);
    v.fill_sp_random(1,1,p);
    Pair<bool> * vpairs;
    int64_t numpair;
    // extract local nonzero data from vector, which in this case is all data
    v.get_local_pairs(&numpair, &vpairs, true);
    // aggregate subcolumns of pairs into bitmasks
    int64_t j = 0;
    while (j < numpair){
      bitmask mask = 0;
      int64_t row = vpairs[j].k / len_bm;
      // accumulate to mask so long as in same subcolumn
      do {
        mask = mask | ((bitmask)1)<<(vpairs[j].k%len_bm);
        j++;
      } while (j < numpair && vpairs[j].k / len_bm == row);
      pairs.push_back(Pair<bitmask>(row+i*mm,mask));
    }
    delete [] vpairs;
  }
  // write data to CTF sparse matrix bulk synchronously from all processors
  A.write(pairs.size(),pairs.data());
 
  return A;
}

/**
 * \brief Given r-by-n matrix A where every element is a bitmask len_bm,
 *        A implicitly stores bool matrix G of size r*len_bm-by-n,
 *        compute B_ij = B_ij + sum_k G_ik G_jk
 *        and     C_ij = C_ij + (sum_k G_ik) + (sum_k G_jk)
 * \param[in] A sparse CTF bitmask matrix, defined as created by generate_random_A() above
 * \param[in,out] B dense CTF matrix of uint64_ts to accmulate to
 * \param[in,out] C dense CTF matrix of uint64_ts to accmulate to
 */
template <typename bitmask>
void jaccard_acc(Matrix<bitmask> & A, Matrix<uint64_t> & B, Matrix<uint64_t> & C){
  B["ij"] += Function<bitmask,bitmask,uint64_t>([](bitmask a, bitmask b){ return (uint64_t)popcount<bitmask>(a&b); })(A["ki"],A["kj"]);

  Vector<uint64_t> v(A.ncol, *A.wrld);
  v["i"] += Function<bitmask,uint64_t>([](bitmask a){ return (uint64_t)popcount<bitmask>(a); })(A["ki"]);
  C["ij"] += v["i"] + v["j"];
}  

/**
 * \brief Compute Jaccard n-by-n dense CTF similarity matrix S,
 *        where S_{ij} = A_{ij}/B_{ij} or 0 if B_{ij} = 0 and A, B are defined jaccard_acc()
 * \param[in] m number of rows in overall bool matrix
 * \param[in] n number of columns in matrix
 * \param[in] p probability of any element being nonzero (even if =1, matrix may have some zeros)
 * \param[in] nbatch number of batches to subdivide the m row into
 * \return S matrix, defined as above
 */
template <typename bitmask>
Matrix<> jaccard_calc_random(int64_t m, int64_t n, double p, int64_t nbatch, World & dw){
  Matrix<> S(n, n, dw);
  Matrix<uint64_t> B(n, n, dw);
  Matrix<uint64_t> C(n, n, dw);

  for (int64_t i=0; i<nbatch; i++){
    int64_t ib = (m/nbatch) + ((m % nbatch) < i);
    Matrix<bitmask> A = generate_random_A<bitmask>(ib, n, p, dw);
    jaccard_acc(A, B, C);
  }
  // subtract intersection from union to get or
  C["ij"] -= B["ij"];
  S["ij"] += Function<uint64_t,uint64_t,double>([](bitmask a, bitmask b){ if (b==0){ assert(a==0); return 0.; } else return (double)a/(double)b; })(B["ij"],C["ij"]);
  
  return S;
}

/**
 * \brief check to make sure each element of S is between 0 and 1
 */
bool is_bounded(Matrix<> & S){
  int64_t num_unbounded = CTF::Function<double,int64_t>([](double s){ return (s>1.) || (s<0.); })(S["ij"]);
  return num_unbounded == 0;
}

/**
 * \brief Test Jaccard similarity matrix computation, by inferring that elements are bounded and that result is the same for bitmasks with 32 and 64 bits
 * \param[in] m number of rows in overall bool matrix
 * \param[in] n number of columns in matrix
 * \param[in] p probability of any element being nonzero (even if =1, matrix may have some zeros)
 * \param[in] nbatch number of batches to subdivide the m row into
 * \return bool true if all tests passed
 */
bool test_jaccard_calc_random(int64_t m, int64_t n, double p, int64_t nbatch){
  World dw(MPI_COMM_WORLD);
  CTF_int::init_rng(dw.rank);
  Matrix<> S32 = jaccard_calc_random<uint32_t>(m, n, p, nbatch, dw);
  CTF_int::init_rng(dw.rank);
  Matrix<> S64 = jaccard_calc_random<uint64_t>(m, n, p, nbatch, dw);

  bool is_bounded_S32 = is_bounded(S32);
  if (!is_bounded_S32){
    if (dw.rank == 0){
      printf("ERROR: uint32_t type calculation of jaccard_calc_random yielded invalid similarities\n");
      return false;
    }
  }
  bool is_bounded_S64 = is_bounded(S64);
  if (!is_bounded_S64){
    if (dw.rank == 0){
      printf("ERROR: uint64_t type calculation of jaccard_calc_random yielded invalid similarities\n");
      return false;
    }
  }

  Matrix<> E(n,n,dw);
  E["ij"] = S32["ij"] - S64["ij"];
  double err32_64 = E.norm2();
  if (err32_64 > 1.e-7*n*n){
    if (dw.rank == 0){
      printf("Similarity matrices disagree, error is %lf\n", err32_64);
    }
    return false;
  } else
    return true;
}

char* getCmdOption(char ** begin,
                   char ** end,
                   const   std::string & option){
  char ** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end){
    return *itr;
  }
  return 0;
}

int main(int argc, char ** argv){
  int rank, np;
  int64_t m, n, nbatch;
  double p;
  int const in_num = argc;
  char ** input_str = argv;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  
  if (getCmdOption(input_str, input_str+in_num, "-m")){
    m = atoi(getCmdOption(input_str, input_str+in_num, "-m"));
    if (m < 0) m = 1023;
  } else m = 1023;

  if (getCmdOption(input_str, input_str+in_num, "-n")){
    n = atoi(getCmdOption(input_str, input_str+in_num, "-n"));
    if (n < 0) n = 15;
  } else n = 15;

  if (getCmdOption(input_str, input_str+in_num, "-p")){
    p = atof(getCmdOption(input_str, input_str+in_num, "-p"));
    if (p < 0) p = .01;
  } else p = .01;
 
  if (getCmdOption(input_str, input_str+in_num, "-nbatch")){
    nbatch = atoi(getCmdOption(input_str, input_str+in_num, "-nbatch"));
    if (nbatch < 0) nbatch = 5;
  } else nbatch = 5;


  if (rank == 0){
    printf("Testing Jaccard similarity matrix construction with %ld-by-%ld k-mer encoding (A) and %ld-by-%ld similarity matrix (S) with nonzero probability in A being p=%lf and number of batches (sets of rows of A) being %ld\n",m,n,n,n,p,nbatch);
  }
  bool pass = test_jaccard_calc_random(m, n, p, nbatch);
  if (rank == 0){
    if (pass)
      printf("Correctness tests passed.\n");
    else
      printf("Correctness tests FAILED!\n");
  }
  MPI_Finalize();
}

