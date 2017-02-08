#include <Rcpp.h>

namespace BasePairs {
  const int A = 1;
  const int C = 2;
  const int G = 3;
  const int T = 4;
  const int Gap = 5;
  const int N = 6;
  inline int GetLargestBPIndex() {
    return N;
  }
  inline Rcpp::CharacterVector GetBPsAtIndex() {
    return Rcpp::CharacterVector::create("A", "C", "G", "T", "-", "N");
  }
}

bool FileExists(const std::string& path);

bool has_suffix(const std::string &str, const std::string &suffix);

inline size_t BPtoIndex(char bp) {
  switch (bp) {
  case 'A':
  case 'a':
    return BasePairs::A;
    break;
  case 'C':
  case 'c':
    return BasePairs::C;
    break;
  case 'G':
  case 'g':
    return BasePairs::G;
    break;
  case 'T':
  case 't':
    return BasePairs::T;
    break;
  case '-':
    return BasePairs::Gap;
    break;
  case 'N':
  case 'n':
    return BasePairs::N;
    break;
  default:
    std::string msg = "Character was not an A, C, G, T, N or -. Check that the data you are loading "
  "(BAM or Reference Fasta) does not contain weird characters.  Strange Character was: ";
     Rcpp::stop(msg + bp);
  }
}
