#include "bithashTest.h"
#include <iostream>

using namespace::std;

CPPUNIT_TEST_SUITE_REGISTRATION(bithashTest);

////////////////////////////////////////////////////////////
// setUp
////////////////////////////////////////////////////////////
void bithashTest::setUp() {

}

////////////////////////////////////////////////////////////
// add
//
// Add a sequence and make sure its there
////////////////////////////////////////////////////////////
void bithashTest::testAdd() {
  const char* nts = "ACGT";
  string seq("CGCGAAAAAACCACA");
  unsigned int iseq[k];
  for (int i = 0; i < k; i++)
    iseq[i] = strchr(nts, seq[i]) - nts;

  bh.add(bh.binary_kmer(seq));
  bh.add(bh.binary_rckmer(seq));

  CPPUNIT_ASSERT(bh.num_kmers() == 2);
  CPPUNIT_ASSERT(bh.check(iseq));
}

////////////////////////////////////////////////////////////
// check long
//
//
// Check a bunch of kmers in a long sequence, by keeping
// track of the current index and making slight changes to it
////////////////////////////////////////////////////////////
void bithashTest::testCheckLong() {
  const char* nts = "ACGT";
  string seq("CGCGAAAAAACCACACGT");
  unsigned int iseq[k+3];
  for (int i = 0; i < k+3; i++)
    iseq[i] = strchr(nts, seq[i]) - nts;

  string seq1("CGCGAAAAAACCACA");
  string seq2("GCGAAAAAACCACAC");
  string seq3("CGAAAAAACCACACG");
  string seq4("GAAAAAACCACACGT");
  bh.add(bh.binary_kmer(seq1));
  bh.add(bh.binary_kmer(seq2));
  //  bh.add(bh.binary_kmer(seq3));
  bh.add(bh.binary_kmer(seq4));

  long long unsigned current;
  CPPUNIT_ASSERT(bh.check(&iseq[0], current));
  CPPUNIT_ASSERT(bh.check(current, iseq[0], iseq[k]));
  CPPUNIT_ASSERT(!bh.check(current, iseq[1], iseq[k+1]));
  CPPUNIT_ASSERT(bh.check(current, iseq[2], iseq[k+2]));
}
