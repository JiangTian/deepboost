#ifndef __DENSITY_HPP__
#define __DENSITY_HPP__

#include <vector>
#include "common.hpp"

class Feature;

class Density {
private:
  std::vector<std::vector<real> > w;
  std::vector<Feature*> phi;
  std::vector<real> expFactors; // = exp(w_t phi(x) )
  real normalizer; // = sum_{x\in S} (exp(w_t phi(x) ) )
public:
  Density() {}
  Density(const std::vector<std::vector<real> > & w,
	  const std::vector<Feature*> & phi,
	  const Dataset & S); //TODO: testing time
  real eval(const Datapoint & X) const;
  real eval(int i) const; // evaluated on S[i] but is faster
};

#endif
