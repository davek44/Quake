#include "qmer_hash.h"
#include <iostream>

////////////////////////////////////////////////////////////////////////////////
// qmer_hash
//
// Form k-mer hash table given a table size, k-mer size, and max count
////////////////////////////////////////////////////////////////////////////////
qmer_hash::qmer_hash(unsigned long long s, unsigned int _k, unsigned int m = 500) {
     if(s & (s-1) != 0) {
	  cerr << "Please use a hash table size that is a power of 2 for double hashing purposes" << endl;
	  exit(1);
     }
     
     size = s;
     table = new unsigned long long[size];
     clear();
     elements = 0;
     k = _k;
     max = m;

     count_mask = 0;
     for(int i = 0; i < (64 - 2*k); i++) {
	  count_mask <<= 1;
	  count_mask |= 1;
     }
     kmer_mask = ~count_mask;

     res_factor = count_mask / m;

     bin_to_char = "ACGT";
}

////////////////////////////////////////////////////////////////////////////////
// ~qmer_hash
//
// De-allocate table
////////////////////////////////////////////////////////////////////////////////
qmer_hash::~qmer_hash() {
     delete[] table;
}


////////////////////////////////////////////////////////////////////////////////
// add
//
// Add q to the kmer's hash table entry.  Will not add if table is full.
////////////////////////////////////////////////////////////////////////////////
void qmer_hash::add(unsigned long long kmer, float q) {
     unsigned long long h1 = dbj_hash(kmer);
     unsigned long long h2 = fnv_hash(kmer);
     if(h2 % 2 == 0)  // force to be odd (Cormen et al. Algorithms)
	  h2 = (h2 + 1) % size;

     bool added = false;
     for(unsigned long long i = 0; i < size && !added; i++) {
	  unsigned long long h = (h1 + i*h2) % size;
	  if(table[h] == 0) {
	       init_entry(h, kmer, q);
	       added = true;

	  } else if(compare_entry(kmer, table[h])) {
	       update_entry(h, q);
	       added = true;
	  }
     }
     
     if(!added)
	  cerr << "qmer_hash is full! element not added!" << endl;
     else
	  elements++;
}


////////////////////////////////////////////////////////////////////////////////
// compare_entry
//
// Compare a kmer to a hash table entry.
////////////////////////////////////////////////////////////////////////////////
bool qmer_hash::compare_entry(unsigned long long kmer, unsigned long long entry) {
     // get kmer from entry
     unsigned long long entry_kmer = entry >> (64 - 2*k);
     return (entry_kmer == kmer);
}


////////////////////////////////////////////////////////////////////////////////
// init_entry
//
// Initialize the kmer and value to the entry h.
////////////////////////////////////////////////////////////////////////////////
void qmer_hash::init_entry(unsigned long long h, unsigned long long kmer, float q) {
     unsigned long long count = (unsigned long long)(res_factor * q);
     if(count > count_mask)
	  count = count_mask;
     table[h] = (kmer << (64 - 2*k)) | count;
     //cout << "Initializing " << h << " to\t" << ull_to_binstr(table[h]) << endl;
}


////////////////////////////////////////////////////////////////////////////////
// update_entry
//
// Add the value q to the entry at h
////////////////////////////////////////////////////////////////////////////////
void qmer_hash::update_entry(unsigned long long h, float q) {
     unsigned long long lkmer = table[h] & kmer_mask;
     unsigned long long count = table[h] & count_mask;

     count += (unsigned long long)(res_factor * q);
     if(count > count_mask)
	  count = count_mask;    

     table[h] = lkmer | count;
}


////////////////////////////////////////////////////////////////////////////////
// print
//
// Print hash table
////////////////////////////////////////////////////////////////////////////////
void qmer_hash::print() {
     unsigned long long count;
     float count_f;
     unsigned long long kmer_ull;
     char* kmer_char = new char[k+1];
     kmer_char[k] = '\0';

     for(unsigned long long i = 0; i < size; i++) {
	  if(table[i] != 0) {
	       //cout << "Entry " << i << ": " << ull_to_binstr(table[i]) << endl;

	       count = table[i] & count_mask;
	       count_f = (float)count / (float)res_factor;
	       
	       kmer_ull = table[i] >> (64 - 2*k);  //so k-mer sits at right	  
	       ull_to_char(kmer_char, kmer_ull);
	       
	       printf("%s\t%.3f\n",kmer_char,count_f);
	  }
     }

     delete[] kmer_char;
}


////////////////////////////////////////////////////////////////////////////////
// ull_to_char
//
// Convert a k-mer as an unsigned long long to a char*
////////////////////////////////////////////////////////////////////////////////
void qmer_hash::ull_to_char(char* kmer_char, unsigned long long kmer_ull) {
     unsigned long long nt;
     for(int i = 0; i < k; i++) {
	  nt = kmer_ull & 3ULL;
	  kmer_char[k-1-i] = bin_to_char[nt];
	  kmer_ull >>= 2;
     }
}


////////////////////////////////////////////////////////////////////////////////
// clear
//
// Clear the table
////////////////////////////////////////////////////////////////////////////////
void qmer_hash::clear() {
     memset(table, 0, size);
     elements = 0;
}


////////////////////////////////////////////////////////////////////////////////
// load
//
// Return load factor of table.
////////////////////////////////////////////////////////////////////////////////
float qmer_hash::load() {
     return (float)elements/float(size);
}


////////////////////////////////////////////////////////////////////////////////
// oat_hash
//
// One-at-a-Time hash
////////////////////////////////////////////////////////////////////////////////
unsigned long long qmer_hash::oat_hash(unsigned long long kmer) {
     unsigned long long nt;
     unsigned long long h = 0;

     for (int i = 0; i < k; i++ ) {
	  nt = kmer & 3;
	  kmer >>=2;

	  h += nt;
	  h += ( h << 10 );
	  h ^= ( h >> 6 );
     }

     h += ( h << 3 );
     h ^= ( h >> 11 );
     h += ( h << 15 );

   return h % size;
}

////////////////////////////////////////////////////////////////////////////////
// fnv_hash
//
// Fowler/Noll/Vo hash
////////////////////////////////////////////////////////////////////////////////
unsigned long long qmer_hash::fnv_hash(unsigned long long kmer) {
     unsigned long long nt;
     unsigned long long h = 14695981039346656037ULL;

     for (int i = 0; i < k; i++ ) {
	  nt = kmer & 3;
	  kmer >>= 2;

	  h = ( h * 1099511628211ULL ) ^ nt;
     }
     
     return h % size;
}

////////////////////////////////////////////////////////////////////////////////
// sax_hash
//
// Shift-Add-XOR hash
////////////////////////////////////////////////////////////////////////////////
unsigned long long qmer_hash::sax_hash(unsigned long long kmer) {
     unsigned long long nt;
     unsigned long long h = 0;

     for (int i = 0; i < k; i++ ) {
          nt = kmer & 3;
          kmer >>= 2;

          h ^= (h << 5) + (h >> 2) + nt;
     }
     
     return h % size;
}

////////////////////////////////////////////////////////////////////////////////
// dbj_hash
//
// Modified Bernstein hash
////////////////////////////////////////////////////////////////////////////////
unsigned long long qmer_hash::dbj_hash(unsigned long long kmer) {
     unsigned long long nt;
     unsigned long long h = 0;

     for (int i = 0; i < k; i++ ) {
          nt = kmer & 3;
          kmer >>= 2;

          h = 67 * h ^ nt;
     }
     
     return h % size;
}

////////////////////////////////////////////////////////////////////////////////
// ull_to_binstr
//
// Convert unsigned long long hash entry to a string of 1's and 0's for
// debugging purposes.
////////////////////////////////////////////////////////////////////////////////
string qmer_hash::ull_to_binstr(unsigned long long ull) {
     string rev;

     for(int i = 0; i < 64; i++) {
	  if(ull & 1 == 1)
	       rev.push_back('1');
	  else
	       rev.push_back('0');

	  ull >>= 1;
     }

     string binstr;
     for(int i = 0; i < 64; i++) {
	  if(i == 2*k)
	       binstr.push_back(',');
	  binstr.push_back(rev[63-i]);
     }

     return binstr;
}
