#include <iostream>
#include <vector>
#include <limits>
#include <cstdio>
#include <cmath>
#include <cstdlib>
using namespace std;

typedef float real;
typedef vector<real> Datapoint;
typedef vector<Datapoint> Dataset;

class Feature; //such as a tree, TODO

struct Density {
  vector<vector<real> > w;
  vector<Feature> phi;
  Density() {}
  Density(const vector<vector<real> > & w, const vector<Feature> & phi)
    :w(w), phi(phi) {}
};

inline real sgn (real x) {
  return (x >= 0.) ? 1. : -1.;
}

real Rm (const Feature & phi) // TODO
real EphiPW(const Feature & phi, const Dataset & S, const Density & pw, int j); // TODO
real EphiS (const Feature & phi, const Dataset & S, int j); // TODO
real Step(const int best_k, const int best_j, const vector<vector<real> > & w, 
const vector<Feature> & phi, const int Lambda, const Dataset & S, const Density & pw, const real beta_k) {
  const int k = best_k;
  const int j = best_j;
  const real wkj = w[k][j];
  const real EphiPWkj = EphiPW(phi[k], S, pw, j);
  const real EphiSkj = EphiS(phi[k], S, j);
  const real pbtp = EphiPWkj + Lambda;
  const real pbtm = EphiPWkj - Lambda;
  const real pbp = EphiSkj + Lambda;
  const real pbm = EphiSkj - Lambda;
  const real e2wL = exp(2. * wkj * Lambda);
  beta = (pbtp * pbm * e2wL - pbp * pbtm) / (pbtp * e2wL - pbtm);
  if (abs(beta) <= beta_k)
    eta = -wkj;
  elseif (beta > beta_k)
    eta = 0.5 / Lambda * log((pbtm * (beta_k - pbp)) / (pbtp * (beta_k - pbm)));
  else
    eta = 0.5 / Lambda * log((pbtm * (beta_k + pbp)) / (pbtp * (beta_k + pbm)));
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
Density DeepMaxent(const Dataset & S, int T) {
  int p = 0;
  vector<int> N;
  vector<vector<real> > w;
  vector<Feature> phi;
  const real tolerance = 1e-3;
  real beta = 0.1; // TODO why ???
  real Lambda = 100.; // TODO
  Density pw;
  
  for (int t = 0; t < T; ++t) {
    real best_d = numeric_limits<real>::min();
    int best_k, best_j = 0;
    for (int k = 0; k < phi.size(); ++k) {
      float beta_k = 2. * Rm(phi[k]) + beta;
      for (int j = 0; j < N[k]; ++j) {
	real d;
	const real wkj = w[k][j];
	const real eps = EphiPW(phi[k], S, pw, j) - EphiS(phi[k], S, j);//TODO optimize
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
	}
      }
    }
    eta = Step(best_k, best_j, w, phi, Lambda, S, pw, beta_k);
    
    w[k][j] += eta;
    
    pw = Density(w, phi);
  }

  return pw
}
