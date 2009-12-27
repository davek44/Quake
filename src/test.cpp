#include <iostream>

using namespace::std;

int main() {
  const char* nts = "ACGT";
  char* i = strchr(nts, 'C');
  if(i) {
    cout << i - nts << endl;
  } else {
    cout << "null" << endl;
  }

  return 0;
}
