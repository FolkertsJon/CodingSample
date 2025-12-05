

// Header file for program that does quadratic interpolation for the flight paths.
#ifndef quadFit_H
#define quadFit_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <math.h>
#include "kUnits.hh"
#include "jonMath.hh"


namespace User
{
  long long int number_of_lines (std::string theFile);
  double* johnTheInterpolator(std::string theFile);
  std::vector<std::vector<double>> ellipticalRadiusTimeTest(double closest, double furthest);
  std::vector<std::vector<double>> ellipticalRadiusTime(double closest, double furthest, double timeStep);
  double ellipticalRadius(double closest, double furthest, double time);
  double* elliptical (double closest, double furthest);
  // This finds the position (distance from the sun) from quadratic constants and time variables
  double positionDays(double time, double* quadConst, double maxTime);
  // hours version
  double positionHours(double time, double* quadConst, double maxTime);
  // minutes version
  double positionMinutes(double time, double* quadConst, double maxTime); 
  // deravative minutes version
  double derivativeMinutes(double time, double* quadConst, double maxTime);

  std::vector<std::pair<double, double>> radiusTimePairsFromFile(std::string filename);
  double interpPosition(std::vector<std::pair<double, double>> pairs, double timeToInterpolate);
}

#endif
