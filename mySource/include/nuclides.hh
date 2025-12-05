

// Header for nuclide processor
#ifndef nuclides_H
#define nuclides_H

#include <string>
#include <iostream>
#include <unordered_map>
#include "kUnits.hh"


namespace User{
  std::string getElementSymbol(int Z);
  
  int getAtomicNumber(const std::string& symbol);
  int getAtomicNumber(int Z);
  
  double getAtomicMass(const std::string& symbol);
  double getAtomicMass(int Z);

  std::pair<int, int> ZAFromPDGCode(const std::string& pdgId);
}

#endif
