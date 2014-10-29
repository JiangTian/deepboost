#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "common.hpp"
#include "feature.hpp"
using namespace std;

//TODO: this could be stored
vector<real> Feature::EphiS(const Dataset & S) const {
  vector<real> ES(this->size(), 0.);
  for (int i = 0; i < S.size(); ++i) {
    vector<real> phix = this->eval(S[i]);
    for (int j = 0; j < this->size(); ++j)
      ES[j] += phix[j];
  }
  for (int j = 0; j < this->size(); ++j)
    ES[j] /= S.size();
  return ES;
}

vector<real> Feature::EphiPW(const Dataset & S, const Density & pw) const {
  vector<real> E(this->size(), 0.);
  for (int i = 0; i < S.size(); ++i) {
    real pwsi = pw.eval(i);
    vector<real> phix = this->eval(S[i]);
    for (int j = 0; j < this->size(); ++j)
      E[j] += pwsi * phix[j];
  }
  return E;
}

// FeatureRaw
vector<real> FeatureRaw::eval(const Datapoint & X) const {
  return vector<real>(1, X[i]);
}

int FeatureRaw::size() const {return 1;}

real FeatureRaw::RademacherComplexity() const {
  return 0;
}

#if 0
// FeatureCategory
vector<real> FeatureCategory::eval(const Datapoint & X) const {
  return (real)(int)(X[i] == category);
}

int FeatureCategory::size() const {return 1;}

real FeatureCategory::RademacherComplexity() const {
  return 0;
}
#endif

// FeatureMonomial2
vector<real> FeatureMonomial2::eval(const Datapoint & X) const {
  return vector<real>(1, X[i] * X[j]);
}

int FeatureMonomial2::size() const {return 1;}

real FeatureMonomial2::RademacherComplexity() const {
  return 0;
}

// FeatureThreshold
vector<real> FeatureThreshold::eval(const Datapoint & X) const {
  return vector<real>(1, (real)(int)(X[i] > threshold));
}

int FeatureThreshold::size() const {return 1;}

real FeatureThreshold::RademacherComplexity() const {
  return 0;
}

// FeatureHinge
vector<real> FeatureHinge::eval(const Datapoint & X) const {
  real out = max((real)0., min((real)1., (X[i] - threshold) / b));
  //real out = (real)(int)(X[i] > threshold) * min((real)1., (X[i] - threshold)/b);
  return vector<real>(1, out);
}

int FeatureHinge::size() const {return 1;}

real FeatureHinge::RademacherComplexity() const {
  return 0;
}
