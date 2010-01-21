#include "bithash.h"
#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace::std;

bithash::bithash() {
  mask = (int)pow(4.0,k) - 1;
}

bithash::~bithash() {
}

////////////////////////////////////////////////////////////
// add
//
// Add a single sequence to the bitmap
////////////////////////////////////////////////////////////
void bithash::add(unsigned long long kmer) {
  bits.set(kmer);
}


////////////////////////////////////////////////////////////
// check
//
// Check for the presence of a sequence in the tree
////////////////////////////////////////////////////////////
bool bithash::check(unsigned kmer[k]) {
  unsigned long long kmermap = 0;
  for(int i = 0; i < k; i++) {
    if(kmer[i] < 4) {
      kmermap <<= 2;
      kmermap |= kmer[i];
    } else
      return false;
  }

  return bits[kmermap];
}

////////////////////////////////////////////////////////////
// check
//
// Check for the presence of a sequence in the tree.
// Pass the kmer map value back by reference to be re-used
////////////////////////////////////////////////////////////
bool bithash::check(unsigned kmer[k], unsigned long long & kmermap) {
  kmermap = 0;
  for(int i = 0; i < k; i++) {
    if(kmer[i] < 4) {
      kmermap <<= 2;
      kmermap |= kmer[i];
    } else
      return false;
  }

  return bits[kmermap];
}

////////////////////////////////////////////////////////////
// check
//
// Check for the presence of a sequence in the tree.
// Pass the kmer map value back by reference to be re-used
////////////////////////////////////////////////////////////
bool bithash::check(unsigned long long & kmermap, unsigned last, unsigned next) {
  if(next >= 4)
    return false;
  else {
    kmermap <<= 2;
    kmermap &= mask;
    kmermap |= next;
  }

  return bits[kmermap];
}


////////////////////////////////////////////////////////////
// file_load
//
// Make a prefix_tree from kmers in the file given that
// occur >= "boundary" times
////////////////////////////////////////////////////////////
void bithash::meryl_file_load(const char* merf, const int boundary) {
  ifstream mer_in(merf);
  string line;
  int count;
  bool add_kmer = false;

  while(getline(mer_in, line)) {
    if(line[0] == '>') {
      // get count
      count = atoi(line.substr(1).c_str());
      //cout << count << endl;
      
      // compare to boundary
      if(count >= boundary) {
	add_kmer = true;
      } else {
	add_kmer = false;
      }

    } else if(add_kmer) {
      // add to tree
      add(binary_kmer(line));

      // add reverse to tree
      add(binary_rckmer(line));
    }
  }
}

////////////////////////////////////////////////////////////
// file_load
//
// Make a prefix_tree from kmers in the file given that
// occur >= "boundary" times
////////////////////////////////////////////////////////////
void bithash::tab_file_load(const char* merf, const int boundary) {
  ifstream mer_in(merf);
  string line;
  int count;

  while(getline(mer_in, line)) {
    if(line[k] != '\t')
      cout << "Kmers are not of expected length " << k << endl;

    // get count
    count = atoi(line.substr(k+1).c_str());
    //cout << count << endl;
      
    // compare to boundary
    if(count >= boundary) {
      // add to tree
      add(binary_kmer(line.substr(0,k)));

      // add reverse to tree
      add(binary_rckmer(line.substr(0,k)));
    }
  }
}

//  Convert string  s  to its binary equivalent in  mer .
unsigned long long  bithash::binary_kmer(const string & s) {
  int  i;
  unsigned long long mer = 0;
  for  (i = 0; i < s.length(); i++) {
    mer <<= 2;
    mer |= binary_nt(s[i]);
  }
  return mer;
}

//  Convert string s to its binary equivalent in mer .
unsigned long long  bithash::binary_rckmer(const string & s) {
  int  i;
  unsigned long long mer = 0;
  for  (i = s.length()-1; i >= 0; i--) {
    mer <<= 2;
    mer |= 3 - binary_nt(s[i]);
  }
  return mer;
}

//  Return the binary equivalent of  ch .
unsigned bithash::binary_nt(char ch) {
  switch  (tolower (ch)) {
  case  'a' : return  0;
  case  'c' : return  1;
  case  'g' : return  2;
  case  't' : return  3;
  }
}


int bithash::num_kmers() {
  return (int)bits.count();
}
