

// Header file for program that does quadratic interpolation for the flight paths.
#ifndef JONMATH_H
#define JONMATH_H

#include<vector>
#include<cmath>
#include<iostream>
#include"kUnits.hh"

namespace User
{

  
  double trapezoidIntegral(double xMin, double yMin, double xMax, double yMax, double intMin, double intMax);
  
  double trapezoidIntegral(double xMin, double yMin, double xMax, double yMax);

  double trapezoidIntegral(std::vector<double> x, std::vector<double> y);

  double trapezoidIntegral(std::vector<double> x, std::vector<double> y, double xMax, double yMax);

  std::pair<long double, long double> wilsonScoreInterval(size_t count, size_t total, long double z);

  double solveKepler(double M, double e);

  double computePositionForKepler(double a, double e, double E);
  
  double computePositionForKeplerChatGPT(double a, double e, double E); // derives r wrongly? No. sqrt(x^2 + y^2) is numericalyl identical to 15 digits of precision in a 100 to 10 solar radii orbit.

  double betheBlochLoss(double z, double a, double beta, double density);
  inline double betheBlochLoss(double z, double beta, double density){
    return  betheBlochLoss(z, 2*z, beta, density); // appriximates a=2z
  }
  
}

#endif
