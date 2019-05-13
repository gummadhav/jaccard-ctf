#include <ctf.hpp>

using namespace CTF;

template <typename bitmask>
uint32_t popcount(bitmask x){
  uint32_t p = 0;
  for (int i=0; i<sizeof(bitmask)*8; i++){
    p += (x>>i) & 1;
  }
}

template <typename bitmask>
Matrix<bitmask> generate_random_A(int64_t m, int64_t n, double p, World & dw){
  Matrix<bitmask> A(m,n,SP,dw);
  std::vector<Pair<bitmask>> pairs((int64_t)((m*p*1.1)*(((double)n)/dw.np +1)));
  World selfw(MPI_COMM_SELF);
  for (int64_t i=dw.rank; i<n; i+=dw.np){
    Vector<bool> v(m,SP,selfw);
    v.fill_sp_random(1,1,p);
    int64_t * inds;
    bool * vals;
    int64_t numpair;
    v.get_local_data(&numpair, &inds, &vals, true);
    int64_t j = 0;
    while (j < numpair){
      bitmask mask = 0;
      int64_t row = inds[j] / (sizeof(bitmask)*8);
      do {
        mask = mask | (1<<(inds[j]%(sizeof(bitmask)*8)));
        j++;
      } while (j < numpair && inds[j] / (sizeof(bitmask)*8) == row);
      printf("inds[%ld]=%ld,mask = %lu\n",j,inds[j],(uint64_t)mask);
      pairs.push_back(Pair<bitmask>(row+i*m,mask));
    }
  }
  A.write(pairs.size(),pairs.data());
 
  return A;
}

template <typename bitmask>
void jaccard_acc(Matrix<bitmask> & A, Matrix<uint64_t> & B, Matrix<uint64_t> & C){
  B["ij"] += Function<bitmask,bitmask,uint64_t>([](bitmask a, bitmask b){ return (uint64_t)popcount(a&b); })(A["ki"],A["kj"]);

  Vector<uint64_t> v(A.ncol, *A.wrld);
  v["i"] += Function<bitmask,uint64_t>([](bitmask a){ return (uint64_t)popcount(a); })(A["ki"]);
  C["ij"] += v["i"] + v["j"];
}  

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
  B.print();
  C.print();
  S["ij"] += Function<uint64_t,uint64_t,double>([](bitmask a, bitmask b){ return (double)a/(double)b; })(B["ik"],C["kj"]);
  
  return S;
}

bool is_bounded(Matrix<> & S){
  S.print();
  int64_t num_unbounded = CTF::Function<double,int64_t>([](double s){ return (s>1.) || (s<0.); })(S["ij"]);
  return num_unbounded == 0;
}

bool test_jaccard_calc_random(int64_t m, int64_t n, double p, int64_t nbatch){
  World dw(MPI_COMM_WORLD);
  Matrix<> S32 = jaccard_calc_random<uint32_t>(m, n, p, nbatch, dw);
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
}

