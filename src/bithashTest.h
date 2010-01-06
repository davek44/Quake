#ifndef BITHASH_TEST_H
#define BITHASH_TEST_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "bithash.h"

class bithashTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(bithashTest);
  CPPUNIT_TEST(testAdd);
  CPPUNIT_TEST(testCheckLong);
  CPPUNIT_TEST_SUITE_END();

 public:
  void setUp();
  void testAdd();
  void testCheckLong();
 private:
  bithash bh;
};

#endif
