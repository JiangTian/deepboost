#ifndef __DENSITY_HPP__
#define __DENSITY_HPP__

#include <vector>
#include "common.hpp"

class Feature;

class Density {
public:
  std::vector<std::vector<real> > w;
  std::vector<Feature*> phi;
  std::vector<real> expFactorsS; // = exp(w_t phi(x) ) on S
  std::vector<real> expFactorsSp; // = exp(w_t phi(x) ) on Sp
  real normalizerS; // = sum_{x\in S} (exp(w_t phi(x) ) ) on S
  real normalizerSp; // = sum_{x\in S} (exp(w_t phi(x) ) ) on Sp
  Dataset Sp; // uniformly drawn at construction
public:
  Density() {}
  Density(const std::vector<std::vector<real> > & w,
	  const std::vector<Feature*> & phi,
	  const Dataset & S, int SpSize); //TODO: testing time
  real evalS(int i) const; // evaluated on S[i] but is faster
  real evalSp(int i) const; // evaluated on Sp[i] but is faster
  void precomputeFactors(const Dataset & S, std::vector<real> & expFactors,
			 real & normalizer);
};

#endif
