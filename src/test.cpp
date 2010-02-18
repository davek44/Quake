#include <iostream>
#include <math.h>

using namespace::std;

int main() {
  int q = 3;
  double p = 1.0-pow(10.0, -q/10.0);
  cout << p << endl;
  int q2 = floor(-10*log(1.0 - p)/log(10));
  cout << q2 << endl;
}
