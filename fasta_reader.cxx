/*
 * This file uses boost library
 * To compile: mpicxx -std=c++0x -O0 -g -Wall -Wno-format -DPMPI -DMPIIO -o fasta_reader fasta_reader.cxx -I/../ctf/include  -L/../ctf/lib -lctf -lblas -L/usr/local/opt/openblas/lib -L/../ctf/scalapack/build/lib -L/usr/local/gfortran/lib -llapack -lblas -lscalapack -lgfortran -lz -lboost_iostreams
 *
 */
#include <ctf.hpp>
#include <string>
#include <cstring>
#include <fstream>
#include <iostream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>


using namespace CTF;

void process_fasta_files(int64_t m, int64_t n, int64_t k, char *infolderPath, char *outfolderPath, const char *listfile, double *perc, World & dw)
{
  // nfiles: number of files this MPI process handles
  int64_t nfiles;
  nfiles = (n / dw.np) + (dw.rank < (n % dw.np));
  int64_t maxfiles;
  // max files are handled by rank 0
  // variable used to sync A.write()s across processes
  maxfiles = (n / dw.np) + (0 < (n % dw.np));
  int64_t maxkmer = *perc * m;

  std::map<char, int> atgc = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};
  FILE *fplist;
  if (listfile != nullptr) {
    fplist = fopen(listfile, "r");
    if (fplist == nullptr && dw.rank == 0) {
      printf("I am unable to open file: %s\n", listfile);
      fflush(stdout);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    char dummy[9000];
    // Each rank starts from its corresponding file
    for (int64_t i = 0; i < dw.rank; i++) {
      fscanf(fplist, "%s", dummy);
    }
  }
  for (int64_t i = 0; i < maxfiles; i++) {
    if (i >= nfiles) {
      continue;
    }

    std::set<int64_t> kmers;
    if (fplist != nullptr) {
      // Read files in lexicographic order
      char dummy[9000];
      fscanf(fplist, "%s", dummy);
      std::string ss = std::string(infolderPath) + "/" + std::string(dummy);
      std::ifstream file(ss.c_str(), std::ios_base::in | std::ios_base::binary);
      if (!file.is_open()) {
        printf("I am rank: %d, I was unable to open file: %s\n", dw.rank, ss.c_str());
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, -1);
      }
      // printf("I am rank: %d, I am opening file: %s\n", dw.rank, ss.c_str());
      boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
      inbuf.push(boost::iostreams::gzip_decompressor());
      inbuf.push(file);
      //Convert streambuf to istream
      std::istream instream(&inbuf);

      ss = std::string(outfolderPath) + "/" + string(dummy).substr(0, string(dummy).size() - 9) + ".txt";
      std::ofstream wfile(ss.c_str());
      if (!wfile.is_open()) {
        printf("I am rank: %d, I was unable to open file: %s\n", dw.rank, ss.c_str());
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, -1);
      }
      // printf("I am rank: %d, I am opening file: %s\n", dw.rank, ss.c_str());

      //Iterate lines
      std::string line;
      std::string pline;
      int64_t ik = 0;
      int64_t kmerv = 0;
      std::map<char, int>::iterator it;
      std::string::size_type il;
      while (std::getline(instream, line)) {
        // std::cout << line << std::endl;
        if (line.find("unitig") != std::string::npos) {
          // header found, process a new read
          ik = 0;
          kmerv = 0;
          continue;
        }
        if (ik != 0) {
          // new line, but not separated by a header, so the read continues
          // prefix the previous line
          line = pline + line;
          il = pline.size();
        }
        else {
          il = 0;
        }
        for(; il < line.size(); il++) {
          it = atgc.find(toupper(line[il]));
          if (it != atgc.end()) {
            kmerv *= 4;
            kmerv += it->second;
            ik++;
            if (ik == k) {
              // store the k-mer value
              // std::cout << kmerv << endl;
              if (kmerv <= maxkmer) {
                // wfile << kmerv << "\n";
                kmers.insert(kmerv);
              }
              ik = 0;
              kmerv = 0;
              il = il - (k - 1);
              assert (il >= 0);
            }
          }
          else {
            // found a non-atgc character
            ik = 0;
            kmerv = 0;
          }
        }
        // store this line as the previous line
        pline = line;
      }
      // write the kmers
      for(auto const& value: kmers) {
        wfile << value << "\n";
      }
      kmers.clear();
      file.close();
      wfile.close();
      for (int64_t i = 0; i < (dw.np - 1); i++) {
        fscanf(fplist, "%s", dummy);
      }
    }
  }
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
  int64_t m, n, k;
  double perc;
  int const in_num = argc;
  char ** input_str = argv;
  char *infolderPath = nullptr;
  char *outfolderPath = nullptr;
  char *listfile = nullptr;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  {
    World dw(MPI_COMM_WORLD);
    
    if (getCmdOption(input_str, input_str+in_num, "-m")){
      m = atoll(getCmdOption(input_str, input_str+in_num, "-m"));
      if (m < 0) m = 1023;
    } else m = 1023;

    if (getCmdOption(input_str, input_str+in_num, "-n")){
      n = atoll(getCmdOption(input_str, input_str+in_num, "-n"));
      if (n < 0) n = 2;
    } else n = 2;
    
    if (getCmdOption(input_str, input_str+in_num, "-k")){
      k = atoll(getCmdOption(input_str, input_str+in_num, "-k"));
      if (k < 0) k = 3;
    } else k = 3;

    if (getCmdOption(input_str, input_str+in_num, "-perc")){
      perc = atof(getCmdOption(input_str, input_str+in_num, "-perc"));
      if (perc < 0) perc = .1;
    } else perc = .1;
 
    if (getCmdOption(input_str, input_str+in_num, "-lfile")) {
       listfile = getCmdOption(input_str, input_str+in_num, "-lfile");
     } else listfile = nullptr;
    
    if (getCmdOption(input_str, input_str+in_num, "-infolderPath")) {
       infolderPath = getCmdOption(input_str, input_str+in_num, "-infolderPath");
     } else infolderPath = nullptr;

    if (getCmdOption(input_str, input_str+in_num, "-outfolderPath")) {
       outfolderPath = getCmdOption(input_str, input_str+in_num, "-outfolderPath");
     } else outfolderPath = nullptr;
    
    if ((listfile == nullptr || infolderPath == nullptr || outfolderPath == nullptr) && rank == 0) {
      printf("Error in the command line parameters");
      fflush(stdout);
      MPI_Abort(MPI_COMM_WORLD, -2);
    }
    process_fasta_files(m, n, k, infolderPath, outfolderPath, listfile, &perc, dw);
  }
  MPI_Finalize();
}

