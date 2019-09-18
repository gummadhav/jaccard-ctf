#include <ctf.hpp>
#include "bitcount.h"
#include "file_reader.h"

using namespace CTF;

/**
 * \brief get number of set bits in uint-like type bitmask
 *        TODO: when possible, can replace this with hardware instructions
 * \param[in] x bitmask cotnaining some 0 and 1 bits
 * \return number of 1 bits in x
 */
/*
template <typename bitmask>
uint32_t popcount(bitmask x){
  uint32_t p = 0;
  for (int i=0; i<sizeof(bitmask)*8; i++){
    p += (x>>i) & 1;
  }
  return p;
}
*/

template <typename bitmask>
uint64_t wfunc(bitmask a, bitmask b) {
  return (uint64_t)popcount<bitmask>(a&b);
}

template <typename bitmask>
Bivar_Function<bitmask,bitmask,uint64_t> * get_jaccard_kernel(){
  Bivar_Kernel<bitmask,bitmask,uint64_t,wfunc> * k = new Bivar_Kernel<bitmask,bitmask,uint64_t,wfunc>();
  k->intersect_only = true;
  return k;
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
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(p, 0.003);
  for (int64_t i=dw.rank; i<n; i+=dw.np){
    // define a bool CTF sparse vector locally, used only to generate a random sparse vector
    Vector<bool> v(m,SP,selfw);
    double p_dist = distribution(generator);
    v.fill_sp_random(1,1,p_dist);
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
  // B["ij"] += Function<bitmask,bitmask,uint64_t>([](bitmask a, bitmask b){ return (uint64_t)popcount<bitmask>(a&b); })(A["ki"],A["kj"]);

  (*get_jaccard_kernel<bitmask>())(A["ki"],A["kj"],B["ij"]);

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
Matrix<> jaccard_calc(int64_t m, int64_t n, double p, int64_t nbatch, World & dw){
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
  Matrix<> S32 = jaccard_calc<uint32_t>(m, n, p, nbatch, dw);
  CTF_int::init_rng(dw.rank);
  Matrix<> S64 = jaccard_calc<uint64_t>(m, n, p, nbatch, dw);

  bool is_bounded_S32 = is_bounded(S32);
  if (!is_bounded_S32){
    if (dw.rank == 0){
      printf("ERROR: uint32_t type calculation of jaccard_calc yielded invalid similarities\n");
      return false;
    }
  }
  bool is_bounded_S64 = is_bounded(S64);
  if (!is_bounded_S64){
    if (dw.rank == 0){
      printf("ERROR: uint64_t type calculation of jaccard_calc yielded invalid similarities\n");
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

template <typename bitmask>
void jacc_calc_from_files(int64_t m, int64_t n, int64_t nbatch, char *gfile, World & dw)
{
    // range from 0 to (2^32 - 1)
    // nbatch should split the range
    // if length is considered to determine nbatch, the length might encompass different kmers from each read (file) depending on whether
    // the kmers are present in the read or not
    //
    // assumption: the kmers in a read/file are sorted
    //
    // each process should read a file/group of files till the specified range
    // squash the zero rows
    // call jacc_acc() to get the result in B and C
    

    // nfiles: number of files this MPI process handles
    int64_t nfiles;
    nfiles = (n / dw.np) + (dw.rank < (n % dw.np));
    int64_t maxfiles;
    // max files are handled by rank 0
    // variable used to sync A.write()s across processes
    maxfiles = (n / dw.np) + (0 < (n % dw.np));
    // maintain file pointers per MPI process
    FILE *fp[nfiles];
    int64_t lastkmer[nfiles];
    for (int64_t i = 0; i < nfiles; i++){
      lastkmer[i] = -1;
    }

    int64_t kmersInBatch = m / nbatch;
    int64_t batchNo = 0;
    int64_t batchStart = batchNo * kmersInBatch;
    int64_t batchEnd = (batchNo + 1) * kmersInBatch - 1;
    // printf("rank: %d nfiles: %lld\n", dw.rank, nfiles);

    Matrix<> S(n, n, dw);
    Matrix<uint64_t> B(n, n, dw);
    Matrix<uint64_t> C(n, n, dw);

    // create matrix m X n
    while (batchNo < nbatch) {
      Timer t_fileRead("File read");
      t_fileRead.start();
      Matrix<int> A(kmersInBatch, n, SP, dw, "hypersparse_A");

      for (int64_t i = 0; i < maxfiles; i++) {
        if (i >= nfiles) {
          A.write(0, nullptr);
          continue;
        }
        // open the file for the first time
        int64_t fileNo = (i * dw.np) + dw.rank;
        if (batchNo == 0) {
          char gfileTemp[10000];
          sprintf(gfileTemp, "%s.%lld.text.annodbg", gfile, fileNo);
          fp[i] = fopen(gfileTemp, "r");
          if (fp[i] == nullptr) {
            printf("I am rank: %d, I was unable to open file: %s", dw.rank, gfileTemp);
            MPI_Abort(MPI_COMM_WORLD, -1);
          }
          // reads to skip the first line
          int64_t dummy;
          fscanf(fp[i], "%lld", &dummy);
          fscanf(fp[i], "%lld", &dummy);
          // printf("rank: %d file opened: %s\n", dw.rank, gfileTemp);
        }

        int64_t nkmersToWrite = 0;
        std::vector<int64_t> gIndex;
        std::vector<int> gData;
        // If there was a last read kmer from the file
        if (lastkmer[i] != -1) {
          gIndex.push_back((lastkmer[i] - batchStart) + fileNo * kmersInBatch);
          gData.push_back(1);
          nkmersToWrite++;
          lastkmer[i] = -1;
        }
        int64_t kmer;
        // read the file till batchEnd or the end of file
        while (fscanf(fp[i], "%lld", &kmer) != EOF) {
          if (kmer > batchEnd) {
            lastkmer[i] = kmer;
            break;
          }
          // write kmer to A
          gIndex.push_back((kmer - batchStart) + fileNo * kmersInBatch);
          gData.push_back(1);
          nkmersToWrite++;
        }
        if (nkmersToWrite != 0) A.write(nkmersToWrite, gIndex.data(), gData.data());
        else A.write(0, nullptr);
      }
      t_fileRead.stop();
      // A.print_matrix();

      Timer t_squashZeroRows("Squash zero rows");
      t_squashZeroRows.start();
      Vector<int> V(n, dw);
      V.fill_sp_random(1, 1, 1);
      Vector<int> R(kmersInBatch, SP, dw);
      // pull out the non-zero rows
      R["i"] = A["ij"] * V["j"];
      int64_t numpair;
      Pair<int> *vpairs;
      R.get_all_pairs(&numpair, &vpairs, true); // R is duplicated across all processes
      // printf("rank: %d numpair: %lld\n", dw.rank, numpair);

      Pair<int> * rowD = new Pair<int>[n];

      int len_bm = sizeof(bitmask) * 8;
      int64_t mm = (numpair + len_bm - 1) / len_bm;
      Matrix<bitmask> J(mm, n, SP, dw, "J");
      // Predefining size to avoid reallocation/copy in push_back(); 
      // if there is imbalance of columns distributed across processes, the corner case should be handled
      Pair<bitmask> *colD = new Pair<bitmask>[mm];
      Pair<int> *colA = new Pair<int>[len_bm];
      // Update J in parallel; the remainder columns are then updated by process 0 alone
      int rem_it = 0;
      uint64_t numColsP = n / dw.np;
      int64_t i = dw.rank; // Column index
      int64_t columnsProcessed = 0;
      while (rem_it < 2) {
        for (; columnsProcessed < numColsP; columnsProcessed++) {
          int64_t rowNo = 0; // To keep track of the new row number
          // Read a column from Matrix A
          int64_t j = 0; // Index into vpairs
          while (j < numpair) {
            int64_t k = 0; // Index into colA; populate mask
            while (j < numpair && k < len_bm) {
              colA[k].k = vpairs[j].k + i * kmersInBatch;
              k++; j++;
            }
            if (rem_it && dw.rank) A.read(0, nullptr); 
            else A.read(k, colA);
            // Store the mask[len_bm] in colD
            bitmask mask = 0;
            for (int64_t l = 0; l < k; l++) {
              // TODO: can read only non-zero data
              if(colA[l].d) {
                mask = mask | ((bitmask)1) << ((colA[l].k % kmersInBatch) % len_bm);
              }
            }
            colD[rowNo].d = mask;
            colD[rowNo].k = rowNo + i * mm;
            rowNo++;
          }
          if (!rem_it) {
            i += dw.np;
            J.write(mm, colD);
          }
          else {
            i++;
            if(dw.rank == 0) J.write(mm, colD);
            else J.write(0, nullptr);
          }
        }
        // handle the remainder columns
        int64_t rem = (n / dw.np) * dw.np;
        i = rem;
        columnsProcessed = rem;
        numColsP = n;
        rem_it++;
      }
      t_squashZeroRows.stop();
      // J.print_matrix();
      // printf("rank: %d J.ncol: %lld J.nrow: %lld\n", dw.rank, J.ncol, J.nrow);
      A.free_self();
      delete [] vpairs;
      Timer t_jaccAcc("jaccard_acc");
      t_jaccAcc.start();
      if (J.ncol != 0 || J.nrow != 0) {
        jaccard_acc(J, B, C);
      }
      t_jaccAcc.stop();
      batchNo++;
      batchStart = batchNo * kmersInBatch;
      batchEnd = (batchNo + 1) * kmersInBatch - 1;
    }
    Timer t_computeS("Compute S");
    t_computeS.start();
    // subtract intersection from union to get or
    C["ij"] -= B["ij"];
    S["ij"] += Function<uint64_t,uint64_t,double>([](bitmask a, bitmask b){ if (b==0){ assert(a==0); return 0.; } else return (double)a/(double)b; })(B["ij"],C["ij"]);
    t_computeS.stop();
    // S.print_matrix();
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
  char *gfile = NULL;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  World dw(MPI_COMM_WORLD);
  
  if (getCmdOption(input_str, input_str+in_num, "-m")){
    m = atoll(getCmdOption(input_str, input_str+in_num, "-m"));
    if (m < 0) m = 1023;
  } else m = 1023;

  if (getCmdOption(input_str, input_str+in_num, "-n")){
    n = atoll(getCmdOption(input_str, input_str+in_num, "-n"));
    if (n < 0) n = 2;
  } else n = 2;

  if (getCmdOption(input_str, input_str+in_num, "-p")){
    p = atof(getCmdOption(input_str, input_str+in_num, "-p"));
    if (p < 0) p = .01;
  } else p = .01;
 
  if (getCmdOption(input_str, input_str+in_num, "-nbatch")){
    nbatch = atoi(getCmdOption(input_str, input_str+in_num, "-nbatch"));
    if (nbatch < 0) nbatch = 1;
  } else nbatch = 1;

  if (getCmdOption(input_str, input_str+in_num, "-f")) {
     gfile = getCmdOption(input_str, input_str+in_num, "-f");
   } else gfile = NULL;
  
  if (gfile != NULL) {
    jacc_calc_from_files<uint32_t>(m, n, nbatch, gfile, dw);
    if (rank == 0) {
      printf("S matrix computed for the specified input dataset\n");
    }
  }

  /*
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
  */
  MPI_Finalize();
}

