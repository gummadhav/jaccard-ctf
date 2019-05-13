#include <ctf.hpp>

using namespace ctf

int popcount(uint64_t x){
  int p = 0;
  for (int i=0; i<64; i++){
    p += (x>>i) & 1;
  }
}

CTF::Matrix generate_random_A(int64_t m, int64_t n, double p, World dw){
  Matrix<uint64_t> A(m,n,SP,dw);
  std::vector<Pair<int64_t,uint64_t>> pairs((int64_t)((m*p)*n/dw.np));
  World selfw(MPI_COMM_SELF);
  for (int64_t i=dw.rank; i<n; i+=dw.np){
    Vector<bool> v(m,SP,selfw);
    v.fill_sp_random(1,1,p);
    int64_t * inds;
    bool * vals;
    int64_t numpair;
    v.get_local_data(&numpair, &inds, &vals, true)
    for (int64_t j=0; j<numpair; j++){
    int64_t j = 0;
    while (j < num_pair){
      uint64_t mask = 0;
      int64_t row = inds[j] / 64;
      do {
        mask += 1<<(inds[j]%64);
        j++;
      } while (j < num_pair && inds[j] / 64 == row);
      pairs.push_back(Pair<uint64_t>(row+i*m,mask));
    }
  }
  A.write(pairs.size(),pairs.data());
 
  return A;
}

void jaccard(CTF::Matrix<uint64_t> A){

  CTF::Matrix<uint64_t> B(A.ncol, A.ncol, *A.wrld)
  CTF::Matrix<uint64_t> C(A.ncol, A.ncol, *A.wrld)
  CTF::Matrix<double> S(A.ncol, A.ncol, *A.wrld)

  B["ij"] += CTF::Function<uint64_t>([](uint64_t a, uint64_t b){ return popcount(a&b); })(A["ik"],A["kj"]);
  
  C["ij"] += CTF::Function<uint64_t>([](uint64_t a, uint64_t b){ return popcount(a|b); })(A["ik"],A["kj"]);
  
  S["ij"] += CTF::Function<uint64_t,uint64_t,double>([](uint64_t a, uint64_t b){ return (double)a/(double)b; }(B["ik"],C["kj"]);
  

}

