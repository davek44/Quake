#include <iostream>
#include <fstream>

using namespace::std;

void read_all(istream & in) {
  string line;
  while(getline(in, line))
    //cout << line << endl;
    int x = 3;
}

int main() {
  streampos orig = cin.tellg();
  cout << "orig: " << orig << endl;
  read_all(cin);
  cout << "after: " << cin.tellg() << endl;
  cin.seekg(orig);
  read_all(cin);

  return 0;
}
