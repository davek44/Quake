#include "ReadTest.h"
#include <iostream>

using namespace::std;

CPPUNIT_TEST_SUITE_REGISTRATION(ReadTest);

////////////////////////////////////////////////////////////
// setUp
////////////////////////////////////////////////////////////
void ReadTest::setUp() {
  int seq[36] = {0};
  string qual(36,'I');
  vector<int> untrusted;
  r = new Read("1", seq, qual, untrusted, 36);
}

////////////////////////////////////////////////////////////
// error_region
//
// Test error_region which computes the region of the read
// that can possibly contain an error given the Read's
// set of untrusted kmers.
////////////////////////////////////////////////////////////
void ReadTest::testRegion() {
  /////////////////////////
  // 1 error tests
  /////////////////////////
  
  // front error - 2
  r->untrusted.clear();
  for(int i = 0; i < 3; i++)
    r->untrusted.push_back(i);
  CPPUNIT_ASSERT(vec_eq(r->error_region(), 0, k-1));

  // back error - 32
  r->untrusted.clear();
  for(int i = 32-k+1; i <= 36-k; i++)
    r->untrusted.push_back(i);
  CPPUNIT_ASSERT(vec_eq(r->error_region(), 36-k, 35));

  /////////////////////////
  // 2 error tests
  /////////////////////////
  
  // 2 front errors - 3, 5
  r->untrusted.clear();
  for(int i = 0; i <= 5; i++)
    r->untrusted.push_back(i);
  CPPUNIT_ASSERT(vec_eq(r->error_region(), 0, k-1));

  // 1 front 1 back - 3, 32
  r->untrusted.clear();
  for(int i = 0; i <= 3; i++)
    r->untrusted.push_back(i);
  for(int i = 32-k+1; i <= 36-k; i++)
    r->untrusted.push_back(i);
  CPPUNIT_ASSERT(vec_eq(r->error_region(), 0, 35));

  // 2 back error - 29, 32
  r->untrusted.clear();
  for(int i = 29-k+1; i <= 36-k; i++)
    r->untrusted.push_back(i);
  CPPUNIT_ASSERT(vec_eq(r->error_region(), 36-k, 35));
}

////////////////////////////////////////////////////////////
// check_trust
//
// Test check_trust to make sure untrusted kmers are
// properly updated
////////////////////////////////////////////////////////////
void ReadTest::testCheckTrust() {
  // make prefix tree
  prefix_tree trusted;
  for(int i = 0; i <= r->read_length-k; i++)
    trusted.add(&r->seq[i]);

  // make change
  r->seq[17] = 1;

  // find untrusted kmers
  r->untrusted.clear();
  for(int i = 0; i <= r->read_length-k; i++) {
    if(!trusted.check(&r->seq[i])) {
      r->untrusted.push_back(i);
    }
  }

  // corrected read
  correction c(17, 0);
  vector<correction*> vc;
  vc.push_back(&c);
  corrected_read cr(vc, r->untrusted, -100, 0);
  
  CPPUNIT_ASSERT(r->check_trust(&cr, &trusted));    
}


////////////////////////////////////////////////////////////
// testRealRead
//
// An H. pylori read that was corrected by Shrec
////////////////////////////////////////////////////////////
bool ReadTest::testRealRead() {
  const char* header = ">@HWI-EAS440_63:6:1:1530:367"
  const char* read = "GTTTTTTTACAGAAGAGAAATTTTTCAGAATGACTA"
  const char* read = "TTCTTGTGATGGAAGAGAAATTTTTCAGAATGACTA";
  

////////////////////////////////////////////////////////////
// vec_eq
//
// Return true if the vector exactly covers the given
// range (inclusively)
////////////////////////////////////////////////////////////
bool ReadTest::vec_eq(vector<int> v, int a, int b) {
  /*
  for(int x = 0; x < v.size(); x++)
    cout << v[x] << endl;
  cout << endl;
  */

  if(v.size() != (b - a + 1))
    return false;

  int i,j;
  bool found;
  for(i = a; i <= b; i++) {
    found = false;
    for(j = 0; j < v.size(); j++) {
      if(i == v[j]) {
	found = true;
	break;
      }
    }
    if(!found)
      return false;
  }
  return true;
}
