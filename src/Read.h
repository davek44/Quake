#ifndef READ_H
#define READ_H

#include "bithash.h"
#include <string>
#include <vector>
#include <fstream>
#include <bitset>

using namespace::std;

//const int bitsize = 22;
const int bitsize = 75;

////////////////////////////////////////////////////////////
// correction
//
// Simple structure for corrections to reads
////////////////////////////////////////////////////////////
class correction {
public:
  correction(short i, short t) {
    index = i;
    to = t;
  };
  // default copy constructor should suffice

  short index;
  //int from;
  short to;
};

////////////////////////////////////////////////////////////
// corrected_read
//
// Simple structure for corrected reads
////////////////////////////////////////////////////////////
class corrected_read {
public:

 corrected_read(vector<correction> & c, bitset<bitsize> & u, float l, short re)
    :untrusted(u) {
    likelihood = l;
    region_edits = re;
    for(int i = 0; i < c.size(); i++)
      corrections.push_back(correction(c[i]));
  };
 corrected_read(bitset<bitsize> & u, float l, short re)
    :untrusted(u) {
    likelihood = l;
    region_edits = re;
  };
  ~corrected_read() {
    /*
    while(corrections.size() > 0) {
      delete corrections.back();
      corrections.pop_back();
    }
    */
  }
  // default destructor should call the correction vector
  // destructor which should call the correction destructor

  vector<correction> corrections; 
  bitset<bitsize> untrusted; // inaccurate until pop'd off queue and processed
  float likelihood;
  short region_edits;
};

////////////////////////////////////////////////////////////
// Read
////////////////////////////////////////////////////////////
class Read {
 public:
  Read(const string & h, const unsigned int* s, const string & q, vector<int> & u, const int read_length);
  ~Read();

  bool trim(int t, ofstream & out);
  bool correct(bithash *trusted, ofstream & out, double (&ntnt_prob)[4][4], bool learning = false);
  bool single_correct(bithash* trusted, ofstream & out, double (&ntnt_prob)[4][4], bool learning = false);
  bool correct_subset(vector<int> untrusted_subset, bithash* trusted, ofstream & out, double (&ntnt_prob)[4][4], bool learning);
  vector<short> error_region(vector<int> untrusted_subset);
  bool check_trust(corrected_read *cr, bithash *trusted);


  string header;
  int read_length;
  int trim_length;
  unsigned int* seq;
  float* prob;
  vector<int> untrusted;
  corrected_read *trusted_read;

  const static float trust_spread_t = .1;
  const static float correct_min_t = .00001;
  const static float learning_min_t = .005;
  static int trim_t;

 private:
  bool untrusted_intersect(vector<int> untrusted_subset, vector<short> & region);
  void untrusted_union(vector<int> untrusted_subset, vector<short> & region);
  void quality_quicksort(vector<short> & indexes, int left, int right);
  string print_seq();
  //string print_corrected(corrected_read* cr);
  string print_corrected(vector<correction> & cor);
  string print_corrected(vector<correction> & cor, int print_nt);

  float likelihood;
};

#endif
