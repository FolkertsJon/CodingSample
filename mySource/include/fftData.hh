// This is the header file for csvParser.cpp
// Jonathan Folkerts Sep 2023
#ifndef fftdata_H
#define fftdata_H


#include<iostream> // system of units
#include <fftw3.h>
#include<kUnits.hh> // system of units
#include<dataStructures.hh> // data structures; inherit vectors and strings from here


namespace User
{
  // takes a data set with time, ch1, ch2, ... and returns frequency, fft1, fft2, ...
  std::vector<std::vector<double>> fftData(std::vector<std::vector<double>> dataSet);

}


#endif
