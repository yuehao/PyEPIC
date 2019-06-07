//
// Created by Yue Hao on 4/15/17.
//
#include <cmath>
#include <string>
#include <vector>
#include <random>

#ifndef EPIC_MATHFUNC_H
#define EPIC_MATHFUNC_H
const double pmass=938272310.0;
const double emass=510999.0;
const double clight=299792458.0;
const double mu0=1.25663706143592e-6;
const double eps0=8.854187817e-12;
const double re0=2.817940325e-15;
const double rp0=1.53469825e-18;
const double sqrtpi=sqrt(M_PI);
const double twopi=2*M_PI;
const std::string tab="\t";

void generateGaussian(std::mt19937& rdgen, std::vector<double>& ranlist,const unsigned long& number, const double& mean, const double& rms);

#endif //EPIC_MATHFUNC_H
