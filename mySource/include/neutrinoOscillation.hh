// This is the header file for neutrinoOscillaion.cpp
// Jonathan Folkerts June 2021
#ifndef variableParser_H
#define variableParser_H

#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <complex>
#include <random>
#include <stdlib.h> 

#include "interp.hh"
#include "dataStructures.hh"
#include "kUnits.hh"
#include "neutrinoOscillation.hh"
#include "fusionModels.hh"

namespace User
{
  Eigen::Matrix3cd hamiltonianStep(User::neutrinoData *theData); // returns the hamiltonian for a particular movemnt through the sun and increments the data; Normal Ordering
  Eigen::Matrix3cd sunSectionHamiltonian ( double minRad, double maxRad , double E);// provides the hamiltonian for a step from minimum radius to masimum radius in the sun at a particular energy; Normal Ordering
  Eigen::Matrix3cd hamiltonianStepIO(User::neutrinoData *theData); // returns the hamiltonian for a particular movemnt through the sun and increments the data; Inverted Ordering
  Eigen::Matrix3cd sunSectionHamiltonianIO( double minRad, double maxRad , double E);// provides the hamiltonian for a step from minimum radius to masimum radius in the sun at a particular energy; Inverted Ordering
  Eigen::Matrix3cd U();//3-d complex double matrix normal ordering
  Eigen::Matrix3cd UIO();//3-d complex double matrix inverse order
  Eigen::Matrix3cd k();//3-d complex double matrix normal order
  Eigen::Matrix3cd kIO();//3-d complex double matrix inverse order
  int factorial(int n);//factorial
  Eigen::Matrix3cd myMatrixExp(Eigen::Matrix3cd theMatrix);// my matrix exponential
  double* solarNeutrinoOscillationProbability(double E);// calculates prob of nu_e,mu, and tau for a particular energy; this oscillates quicklky down to 1e-6 MeV changes at 1 MeV
  double* solarNeutrinoSurvival(double E, double radius);// calculates prob of nu_e,mu, and tau for a particular energy; this oscillates quicklky down to 1e-6 MeV changes at 1 MeV
  
  
  
}








#endif
