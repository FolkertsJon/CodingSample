// This is the header file for fusionModels.cpp
// Jonathan Folkerts June 2021
#ifndef fusionModels_H
#define fusionModels_H

#include <math.h>
#include <random> 
#include <iostream> 
#include <stdlib.h> 

#include "dataStructures.hh"
#include "kUnits.hh"

namespace User
{
  void shellFusionModel(double radius, User::neutrinoData *theData); // puts fusion shell model into neutrino data
  
  
  
}


#endif
