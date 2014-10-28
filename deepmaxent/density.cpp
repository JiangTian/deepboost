#include <iostream>
#include <vector>
#include <cmath>
#include "common.hpp"
#include "density.hpp"
#include "feature.hpp"
using namespace std;

Density::Density(const vector<vector<real> > & w,
		 const vector<Feature*> & phi,
		 const Dataset & S) //TODO: testing time
  :w(w), phi(phi), normalizer(0.), expFactors(S.size(), 0.) {
  for (int i = 0; i < S.size(); ++i) {
    real t = 0.;
    for (int j = 0; j < phi.size(); ++j) {
      vector<real> phix = phi[j]->eval(S[i]);
      for (int k = 0; k < phix.size(); ++k)
	t += w[j][k] * phix[k];
    }
    expFactors[i] = exp(t);
  }
  for (int i = 0; i < expFactors.size(); ++i)
    normalizer += expFactors[i];
}

real Density::eval(const Datapoint & X) const {
  // TODO this can be cached
  real t = 0.;
  for (int j = 0; j < phi.size(); ++j) {
    vector<real> phix = phi[j]->eval(X);
    for (int k = 0; k < phix.size(); ++k)
      t += w[j][k] * phix[k];
  }
  return exp(t) / normalizer;
}

real Density::eval(int i) const {
  cout << i << " " << expFactors.size() << endl;
  return expFactors[i] / normalizer;
}
