#include <iostream>
#include "bithash.h"
#include <string>

using namespace::std;

int main() {
  bithash bh;
  
  string s("AAAACCGGTTGCA");
  bh.add(bh.binary_kmer(s));
}
