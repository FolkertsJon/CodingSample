// Header file for the conversion from root files out of
// the fast monte carlo to a root file containing the
// luminosity of each simulated mission source code
// JF - 2/2024
#ifndef FLIGHTTOLUMINOSITY_H
#define FLIGHTTOLUMINOSITY_H

#include <iomanip> // setprecision
#include <cmath>

#include "TNtuple.h"
#include "TFile.h"

#include "galliumInteraction.hh"
#include "dataStructures.hh"

namespace User
{
  
  
  // Global variables
  extern double nOrbits;
  extern double timeOfFlight;
  extern double closest;
  extern double furthest;
  extern double kgGa;
  extern bool ellipticalOrbit;
  extern std::string targetName;
  extern double nuPerMeVFusion;
  extern double maxRadius;



  int singleFile(std::string filename, std::string location, TNtuple *toBeFilled);// what it says on the tin; outputs to terminal
  
  int directCalc(std::string filename, std::string location, TNtuple *toBeFilled);// what it says on the tin; outputs to terminal; for direct flight calcs

  int manyFiles(std::string filePath, std::string fileNameBase, unsigned int numberOfFiles, TNtuple *toBeFilled);// operates on many files; outputs to ROOT file

  void parseFile(std::string filename, std::string location, TNtuple *toBeFilled);
  
  void parseFileDirectCalc(std::string filename, std::string location, TNtuple *toBeFilled);

  neutrinoTarget targetNameToData(std::string targetName);

  double radiusToTimeScalar(TNtuple* radiusTuple, double maxRadius);

  std::pair<double,double> countNeutrinos(TNtuple* neutrinoTuple, double totalCrossSection);

  std::pair<double, double> neutrinoCountAndUncertainty();

  double nuPerMeVFusionFiller();
  
}

#endif
