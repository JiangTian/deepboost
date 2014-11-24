#include "densityGrid.hpp"
#include <iostream>
#include <cmath>
#include <limits>
#include "common.hpp"
#include "feature.hpp"
using namespace std;

DensityGrid::DensityGrid(const vector<vector<real> > & w,
			 const vector<Feature*> & phi,
			 const Dataset & S, const Dataset & worldGrid)
  :w(w), phi(phi), normalizer(0.) {
  // compute the normalizer
  real largestFactor = -numeric_limits<real>::max();
  for (int i = 0; i < worldGrid.size(); ++i) {
    real t = 0.;
    for (int j = 0; j < phi.size() ++j) {
      vector<real> phix = phi[j]->eval(worldGrid[i]);
      for (int k = 0; k < phix.size(); ++k)
	t += w[j][k] * phix[k];
    }
    
}
