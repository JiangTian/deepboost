#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <cstdlib>
#include "common.hpp"
#include "density.hpp"
#include "feature.hpp"
using namespace std;

Density::Density(const vector<vector<real> > & w,
		 const vector<Feature*> & phi,
		 const Dataset & S, int SpSize) //TODO: testing time
  :w(w), phi(phi), normalizerS(0.), normalizerSp(0.),
   expFactorsS(S.size(), 0.), expFactorsSp(SpSize, 0.), Sp(SpSize) {
  // precompute on S
  precomputeFactors(S, expFactorsS, normalizerS); // we don't need that!

  // generate Sp
  //  get bounds. TODO: get these bounds once and for all
  int inputSize = S[0].size();
  vector<real> xmin(inputSize, numeric_limits<real>::max());
  vector<real> xmax(inputSize, -numeric_limits<real>::max());
  for (int i = 0; i < S.size(); ++i)
    for (int j = 0; j < S[i].size(); ++j) {
      xmin[j] = min(xmin[j], S[i][j]);
      xmax[j] = max(xmax[j], S[i][j]);
    }
  //  sample Sp between bounds
  for (int i = 0; i < SpSize; ++i) {
    Sp[i] = Datapoint(inputSize);
    for (int j = 0; j < inputSize; ++j)
      Sp[i][j] = uniform(xmin[j], xmax[j]);
  }

  // precompute on Sp
  precomputeFactors(Sp, expFactorsSp, normalizerSp);
}

void Density::precomputeFactors(const Dataset & S, vector<real> & expFactors,
				real & normalizer) {
  for (int i = 0; i < S.size(); ++i) {
    real t = 0.;
    for (int j = 0; j < phi.size(); ++j) {
      vector<real> phix = phi[j]->eval(S[i]);
      for (int k = 0; k < phix.size(); ++k)
	t += w[j][k] * phix[k];
    }
    expFactors[i] = exp(t);
    if (expFactors[i] > numeric_limits<real>::max()) {
      cerr << "The factor is too big for machine precision. "
	   << "Likely that w became too large" << endl;
      exit(0);
    }
  }
  for (int i = 0; i < expFactors.size(); ++i)
    normalizer += expFactors[i];
}

real Density::evalS(int i) const {
  return expFactorsS[i] / normalizerS;
}

real Density::evalSp(int i) const {
  return expFactorsSp[i] / normalizerSp;
}
