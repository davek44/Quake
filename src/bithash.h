#ifndef BITHASH_H
#define BITHASH_H

#include <string>
#include <cmath>
#include <bitset>
using namespace::std;

const int k = 15;
const int bitssize = 1073741824;  // i.e. 4^k

class bithash {
 public:
  bithash();
  ~bithash();
  void add(long long unsigned kmer);
  bool check(unsigned kmer[k]);
  bool check(unsigned kmer[k], long long unsigned & kmermap);
  bool check(long long unsigned & kmermap, unsigned last, unsigned next);
  void file_load(const char* merf, const int boundary);
  long long unsigned binary_kmer(const string &s);
  long long unsigned binary_rckmer(const string &s);
  int num_kmers();

 private:
  unsigned binary_nt(char ch);

  bitset<bitssize> bits;
  unsigned long long mask;
};

#endif
