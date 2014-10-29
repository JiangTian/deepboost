#include <iostream>
#include <cmath>
#include <cstdio>
#include <vector>
#include <cstdlib>
#include <cfloat>
#include "deepmaxent.hpp"
using namespace std;

real normal() {
  real x = 0.;
  int n = 10;
  for (int i = 1; i < n; ++i)
    x += 2. * ((real)rand() / (real)RAND_MAX) - 1.;
  return x / sqrt(n);
}

int main(int argc, char* argv[]) {
  Dataset S;
  int nDatapoints = 42;
  int inputDim = 10;

  for (int i = 0; i < nDatapoints; ++i) {
    Datapoint p(inputDim);
    for (int j = 0; j < inputDim; ++j)
      p[j] = normal();
    S.push_back(p);
  }

  Density pw = DeepMaxent(S, 10000);

#if 0
  for (int i = 0; i < pw.w.size(); ++i)
    for (int j = 0; j < pw.w[i].size(); ++j)
      cout << pw.w[i][j] << " ";
  cout << endl;
#endif

  return 0;
}
