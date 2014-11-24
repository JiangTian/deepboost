#ifndef __DENSITY_GRID_HPP__
#define __DENSITY_GRID_HPP__

#include <vector>
#include "common.hpp"

class Feature;

class DensityGrid {
public:
  std::vector<std::vector<real> > w;
  std::vector<Feature*> phi;
  real normalizer;
public:
  DensityGrid() {}
  DensityGrid(const std::vector<std::vector<real> > & w,
	  const std::vector<Feature*> & phi,
	  const Dataset & S, const Dataset & worldGrid);
  real eval(const Datapoint & X) const;
  real evalS(int i) const; // evaluated on S[i] but is faster
  real loss(const Dataset & S, real beta) const;
};

#endif
