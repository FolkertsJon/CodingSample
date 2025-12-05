// This is the header file for oscillations.cpp
// Jonathan Folkerts May 2021
#ifndef oscillations_H
#define oscillations_H


#include<kUnits.hh> // system of units
#include<dataStructures.hh> // data structures
#include<galliumInteraction.hh> // flux approximations
#include<interp.hh>// interpolation


namespace User
{

  std::vector<User::energyProbPair> readData(const std::string& filename);
  
  double nu_e_Earth(const std::vector<User::energyProbPair>& data, double energy);


}

#endif
