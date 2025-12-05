// This is the header file for csvParser.cpp
// Jonathan Folkerts Sep 2023
#ifndef csvparser_H
#define csvparser_H


#include<iostream> 
#include<fstream> 

#include<kUnits.hh> // system of units
#include<dataStructures.hh> // data structures


namespace User
{
  // returns a vector holding 2d data: Time, ch1, ch2, ch3, ch4
  std::vector<std::vector<double>> tektronicsCSVReader(std::string filename);
  std::vector<std::vector<double>> tektronicsCSVReaderDigital(std::string filename, int digitalChannelCount);

  // returns the vertical scales of the csv
  std::vector<double> tektronicsScaleReader(std::string filename);


  // Read TEK files in new format
  oscilloscopeData tekToOscopeData(std::string filename);

  // Smooths the tek data for derivative
  std::vector<double> generateGaussianKernel(int kernelSize, double sigma);
  std::vector<double> tekGausSmooth(oscilloscopeData data, size_t channelNumber, double stDev);

  // Returns the time derivative of the oscilloscope data
  std::vector<double> tekDerivative(oscilloscopeData data, size_t channel);

  // Returns approximate index of peaks from a smoothed derivative
  std::vector<std::vector<size_t>> peakLocationsFromDerivative(User::oscilloscopeData derivativeData, User::oscilloscopeData data);

  bool peakPassesEnergyCutIntegral(size_t peakIdx, double min_in_Vs, double max_in_Vs, User::oscilloscopeData data, size_t channel);


  bool peakPassesEnergyCutVoltage(size_t peakIdx, double min_in_V, double max_in_V, User::oscilloscopeData data, size_t channel);
}

#endif
