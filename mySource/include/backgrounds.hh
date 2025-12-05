

// Header file for the gallium interaction source code
// JF - 2/2021
#ifndef BACKGROUNDS_H
#define BACKGROUNDS_H

#include "kUnits.hh"

namespace User
{
  // This function determines the neutron cross section
  double neutronCrossSection (double radius, double inclination, double energy);

  // This program uses approximations from
  // https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2011SW000732
  // to find the neutron flux. Assuming each flare has approximately the 
  // same overall flux, we can find the flux from one 
  double neutronFlux (double radius, double inclination, double energy);


  // This program takes the detector's position (radius, inclination)
  // and returns a double holding the cross section of a solar wind 
  // background interaction. We want to code softly enough that the 
  // function can run in double-pulse and single pulse modes.
  double neutronTargetNumber(); // 10 times gallium target number for now
  double protonTargetNumber(); // ~assuming the target is 400 pounds of tungsten

  double solarBackground (double radius, double inclination, double detectorDiameter);

  // approximating from https://www.researchgate.net/publication/235341190_SUPERCONDUCTING_TECHNOLOGIES_FOR_THE_ADVANCED_ACCELERATOR_APPLIATIONS_PROGRAM
  // we have that protons below 30 MeV will be stopped by 0.1 cm of tunsten

  // this program returns the GCR flux at an energy from
  // https://alteaspace.wordpress.com/2011/11/27/galactic-cosmic-rays-gcr/
  // which gives approximate GCR flux.
  double gcrProtonFlux (double radius, double inclination, double energy);

  // total flux
  double gcrTotalProtonFlux (double radius, double inclination);

  // This returns the approximate total cross section of the same
  double protonTotalCrossSection (double radius, double inclination);

  // This program returns the rate of cosmic background in units of s^-1
  double cosmicBackground (double radius, double inclination);
    
  // This program returns a probability of a radiological decay within a given timespan
  double radiologicalBackground ();
}





















#endif
