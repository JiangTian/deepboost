#include <iostream>
#include <cmath>
#include <cstdio>
#include <vector>
#include <cstdlib>
#include <random>
#include "deepmaxent.hpp"
using namespace std;

int main(int argc, char* argv[]) {
  Dataset S;
  int nDatapoints = 42;
  int inputDim = 10;

  default_random_engine generator;
  normal_distribution<real> normal(0., 1.);

  for (int i = 0; i < nDatapoints; ++i) {
    Datapoint p(inputDim);
    for (int j = 0; j < inputDim; ++j)
      p[j] = normal(generator);
    S.push_back(p);
  }

  Density pw = DeepMaxent(S, 1);

  return 0;
}
