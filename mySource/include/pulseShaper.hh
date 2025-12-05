// This is the header file for pulseShaper.cpp
// Jonathan Folkerts June 2022
#ifndef pulseShaper_H
#define pulseShaper_H

#include<vector>
#include<algorithm>
#include<iostream>

#include<kUnits.hh>
#include<dataStructures.hh>

namespace User
{
  
  std::vector<double> *vetoPulseShaper( std::vector<double>* pulseIn, std::vector<double>* pulseTime);// This function will take in a vector with the pulse, and output the shaped pulse

  std::vector<double> *vetoRCTimeAndMax( std::vector<double>* pulseIn, std::vector<double>* pulseTime, std::vector<bool>* vetoArr);// This spits out the RC time to go from max/e to the max and the maximum value

std::vector<double> *vetoRCTimeAndMax( std::vector<double>* pulseIn, std::vector<double>* pulseTime);// This spits out the RC time to go from max/e to the max and the maximum value

  
  std::vector<double> *signalMaxTimeAndMax( std::vector<double>* pulseIn, std::vector<double>* pulseTime); // spits out a vector with the maximum value of the input signal and the time of that maximum

  
  std::vector<bool> *produceVeto( std::vector<double>* vetoPulse, double timeStep); // returns a vector with true (veto at this value) or false (don't veto)


  double pulseArea(std::vector<double>* pulseIn, double stepSize); // spits out the area of a pulse from the 1/e time before the max to the 1/e time after the max

}


#endif
