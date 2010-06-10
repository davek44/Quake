#include  <string>
#include  <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
// reduce-kmers
//
// Accept sorted input from count-kmers and reduce by summing counts with the
// same header (which will always be adjacent).
////////////////////////////////////////////////////////////////////////////////
int  main (int argc, char * argv [])
{
  string mykmer, kmer;
  int mycount, count;
  
  // read initial line
  cin >> mykmer;
  cin >> mycount;
  
  while(cin >> kmer) {
    cin >> count;
    
    // if same as last, increment count
    if(kmer == mykmer) {
      mycount += count;
      
      // else print last, initialize new
    } else {
      cout << mykmer << "\t" << mycount << endl;
      mykmer = kmer;
      mycount = count;
    }
  }
  
  // print last
  cout << mykmer << "\t" << mycount << endl;
  
  return 0;
}
