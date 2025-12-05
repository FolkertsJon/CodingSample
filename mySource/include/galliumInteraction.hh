

// Header file for the gallium interaction source code
// JF - 2/2021
#ifndef GALLIUMINTERACTION_H
#define GALLIUMINTERACTION_H

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <omp.h>
//#include <cmath>
#include <chrono>
#include <thread>

#include "TRandom3.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TROOT.h"
#include "TGraph.h"

#include "kUnits.hh"
#include "jonMath.hh"
#include "dataStructures.hh"
#include "oscillations.hh"

namespace User
{
  // NEUTRINOS BY PROCESS AND ENERGY


  extern bool excitedOnly;
  extern bool doOscillations;
  extern std::vector<User::energyProbPair> oscillationData;

  
  inline double galliumThreshold(){
    if(excitedOnly) return 408*kUnits::keV; // threshold for Ga71 with x-ray
    return 233.2*kUnits::keV; // bare threshold
  }
  
  // NEUTRINOS BY PROCESS AND ENERGY
  double powerIntegral(double xMin, double yMin, double xMax, double yMax, double intMin, double intMax);
  //double trapezoidIntegral(double xMin, double yMin, double xMax, double yMax, double intMin, double intMax);

  TGraph *oscillationProbabilityGraph();
  extern TGraph *globalOscillationGraph;
  double oscillationProbability(double energy);

  double hepFlux (double eMin, double eMax);
  double B8Flux (double eMin, double eMax);
  double F17Flux (double eMin, double eMax);
  double O15Flux (double eMin, double eMax);
  double N13Flux (double eMin, double eMax);
  double ppFlux (double eMin, double eMax);
  double pepFlux(double eMin, double eMax);
  double Be7Flux(double eMin, double eMax);

  extern std::vector<double> *Fluxes;
  std::vector<double>* fluxesFiller();
  std::vector<double>* fluxReturner();
  std::vector<double>* fluxesFiller(double eMin);
  std::vector<double>* fluxReturner(double eMin);
  double totalFlux();
  extern double *totalFLuxPointer;
  

  // Cross sections for processes

  double ppCrossSec(double eMin, double eMax);
  double pepCrossSec(double eMin, double eMax);
  double Be7CrossSec(double eMin, double eMax);
  double N13CrossSec(double eMin, double eMax);
  double O15CrossSec(double eMin, double eMax);
  double F17CrossSec(double eMin, double eMax);
  double B8CrossSec(double eMin, double eMax);
  double hepCrossSec(double eMin, double eMax);
  double bahcallEffCrossSection();

  // Generate TGraphs of flux vs energy
  
  TGraph* hepEnergyGraphMaker();
  TGraph* ppEnergyGraphMaker();
  TGraph* N13EnergyGraphMaker();
  TGraph* O15EnergyGraphMaker();
  TGraph* F17EnergyGraphMaker();
  TGraph* B8EnergyGraphMaker();
  TGraph* pepEnergyGraphMaker();
  TGraph* Be7EnergyGraphMaker();

  // Generate histograms of flux vs energy

  TH1D* hepEnergyHistogramMaker();
  TH1D* ppEnergyHistogramMaker();
  TH1D* N13EnergyHistogramMaker();
  TH1D* O15EnergyHistogramMaker();
  TH1D* F17EnergyHistogramMaker();
  TH1D* B8EnergyHistogramMaker();

  // Generate random Energy for processes

  double ppEnergy(TH1D* hist);
  double pepEnergy();
  double Be7Energy();
  double N13Energy(TH1D* hist);
  double O15Energy(TH1D* hist);
  double F17Energy(TH1D* hist);
  double B8Energy(TH1D* hist);
  double hepEnergy(TH1D* hist);

  // The Rate per atom per second at 1 AU
  double neutrinoRate();
  double neutrinoRate(double eMin, double eMax);
  double neutrinoRateByProcess(double eMin, double eMax);

  // Oscillation Models
  double linOscFactor(double Radius);

  // random energy for neutrino
  double neutrinoEnergy();

  // Check values
  void checkGalliumCrossSection();
  
}

#endif
