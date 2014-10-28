#include <iostream>
#include <vector>
#include <limits>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include "common.hpp"
#include "feature.hpp"
#include "density.hpp"
using namespace std;

inline real sgn (real x) {
  return (x >= 0.) ? 1. : -1.;
}

real Step(const int best_k, const int best_j, const vector<vector<real> > & w, 
	  const vector<Feature*> & phi, const int Lambda, const Dataset & S,
	  const Density & pw, const real beta_k) {
  const int k = best_k;
  const int j = best_j;
  const real wkj = w[k][j];
  const real EphiPWkj = (phi[k]->EphiPW(S, pw))[j];
  const real EphiSkj = (phi[k]->EphiS(S))[j];
  const real pbtp = EphiPWkj + Lambda;
  const real pbtm = EphiPWkj - Lambda;
  const real pbp = EphiSkj + Lambda;
  const real pbm = EphiSkj - Lambda;
  const real e2wL = exp(2. * wkj * Lambda);
  const real beta = (pbtp * pbm * e2wL - pbp * pbtm) / (pbtp * e2wL - pbtm);
  if (abs(beta) <= beta_k)
    return -wkj;
  else if (beta > beta_k)
    return 0.5 / Lambda * log((pbtm * (beta_k - pbp)) / (pbtp * (beta_k - pbm)));
  else
    return 0.5 / Lambda * log((pbtm * (beta_k + pbp)) / (pbtp * (beta_k + pbm)));
}

// S : dataset
// T : number of iterations
// N : array. N[k] is the dimension of the output of tree k
// w : (w_t in paper) : vector of vectors of weights
// d (best_d) : cf paper (d_kj in paper)
// phi : trees
// tolerance : epsilon for comparison with 0
// beta_k : 
// eps : (epsilon_t-1,k,j in paper)
// epsK : array of eps for fixed t and k
// phibar, phibarT : 
// Lambda : total range of the phi (TODO)
// eta : size of the update of w
// inputDim : dimension of the input
Density DeepMaxent(const Dataset & S, int T) {

  cout << "1" << endl;
  
  // initialization
  int p = 0;
  vector<vector<real> > w;
  vector<Feature*> phi;
  const real tolerance = 1e-3;
  real beta = 0.1; // TODO why ???
  real Lambda = 100.; // TODO
  int inputDim = S[0].size();
  // TODO: check that S is consistent (all samples have same dimension)

  cout << "2" << endl;

  // adding all possible features
  //  adding raw features
  for (int i = 0; i < inputDim; ++i) {
    Feature* newFeature = new FeatureRaw(i);
    phi.push_back(newFeature);
  }

  cout << "3" << endl;

  //  adding categorical features
#if 0
  for (int i = 0; i < inputDim; ++i) {
    //set category //TODO!!
    S[i]
    Feature* newFeature = new FeatureCategory(i, category);
    phi.push_back(newFeature);
  }
#endif
  //  adding threshold features
  {
    vector<real> sortedInput;
    for (int i = 0; i < inputDim; ++i) {
      sortedInput.clear();
      for (int j = 0; j < S.size(); ++j)
	sortedInput.push_back(S[j][i]);
      sort(sortedInput.begin(), sortedInput.end());
      for (int j = 1; j < sortedInput.size(); ++j) {
	real threshold = 0.5 * (sortedInput[j-1] + sortedInput[j]);
	Feature* newFeature = new FeatureThreshold(i, threshold);
	// TODO (or not): this is never destroyed
	phi.push_back(newFeature);
      }
    }
  }

  cout << "4" << endl;
  
  //  adding hinge features
  {
    real b = 1;
    vector<real> sortedInput;
    for (int i = 0; i < inputDim; ++i) {
      sortedInput.clear();
      for (int j = 0; j < S.size(); ++j)
        sortedInput.push_back(S[j][i]);
      sort(sortedInput.begin(), sortedInput.end());
      for (int j = 1; j < sortedInput.size(); ++j) {
        real threshold = 0.5 * (sortedInput[j-1] + sortedInput[j]);
        Feature* newFeature = new FeatureHinge(i, threshold, b);
        phi.push_back(newFeature);
      }
    }
  }
  
  cout << "5" << endl;

  //  adding monomial2 features											      
  for (int i = 0; i < inputDim; ++i)
    for (int j = 0; j <= i; ++j) {
      Feature* newFeature = new FeatureMonomial2(i, j);
      phi.push_back(newFeature);
    }

  cout << "6" << endl;

  // initial distribution (all w's = 1)
  for (int i = 0; i < phi.size(); ++i)
    w.push_back(vector<real>(phi[i]->size(), 1.));
  Density pw(w, phi, S);
  
  // main loop (changing the w's)
  for (int t = 0; t < T; ++t) {
    real best_d = numeric_limits<real>::min();
    real best_beta_k = 0.;
    int best_k, best_j = 0;
    cout << " 1" << endl;
    for (int k = 0; k < phi.size(); ++k) {
      cout << " 2" << endl;
      float beta_k = 2. * phi[k]->RademacherComplexity() + beta;
      cout << " 3" << endl;
      const vector<real> EphiPW = phi[k]->EphiPW(S, pw);
      cout << " 4" << endl;
      const vector<real> EphiS = phi[k]->EphiS(S);
      cout << " 5" << endl;
      for (int j = 0; j < phi[k]->size(); ++j) {
	cout << " 6" << endl;
	real d;
	const real wkj = w[k][j];
	const real eps = EphiPW[j] - EphiS[j];
	cout << " 7" << endl;
	if (abs(wkj) > tolerance)
	  d = beta_k * sgn(wkj) + eps;
	else if (abs(eps) <= beta_k)
	  d = 0.;
	else
	  d = - beta_k * sgn(wkj) + eps;
	if (d > best_d) {
	  best_d = d;
	  best_k = k;
	  best_j = j;
	  best_beta_k = beta_k;
	}
      }
    }
    real eta = Step(best_k, best_j, w, phi, Lambda, S, pw, best_beta_k);
    
    w[best_k][best_j] += eta;
    
    pw = Density(w, phi, S);
  }

  cout << "7" << endl;

  return pw;
}
