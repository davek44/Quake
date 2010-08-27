#include <string>

using namespace::std;

////////////////////////////////////////////////////////////////////////////////
// qmer_hash
//
// Hash table for q-mer counts that uses only 64 bits per kmer and count.
//
// Uses open addressing via double hashing.
//
// Based on what k is set to, the remainder of the 64 bits are used to store
// the count up to a max provided by the user.  This is done by multiplying
// every float count by a resolution factor and then storing the result int.
// Then to print the count, you divide back by that resolution factor.
//
// No dynamic re-sizing.  Will fail if it fills.
////////////////////////////////////////////////////////////////////////////////
class qmer_hash {
 public:
  qmer_hash(unsigned long long s, unsigned int _k, unsigned int m);
  ~qmer_hash();
  void add(unsigned long long kmer, float q);
  void print();
  void clear();
  float load();
  
 private:
  unsigned long long fnv_hash(unsigned long long kmer);
  unsigned long long oat_hash(unsigned long long kmer);
  unsigned long long dbj_hash(unsigned long long kmer);
  unsigned long long sax_hash(unsigned long long kmer);
  void init_entry(unsigned long long h, unsigned long long kmer, float q);
  bool compare_entry(unsigned long long kmer, unsigned long long entry);
  void update_entry(unsigned long long h, float q);
  void ull_to_char(char* kmer_char, unsigned long long kmer_ull);
  string ull_to_binstr(unsigned long long ull);  

  unsigned long long * table;
  unsigned long long size;
  unsigned long long elements;
  unsigned int max;
  unsigned int k;
  unsigned int res_factor;

  unsigned long long kmer_mask;
  unsigned long long count_mask;
  const char* bin_to_char;
};
