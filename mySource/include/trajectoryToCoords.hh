


// Header file for the gallium interaction source code
// JF - 2/2021
#ifndef trajectoryToCoords_H
#define trajectoryToCoords_H

#include<fstream>
#include<string>
#include<iostream>
#include<sstream>
#include <vector>
#include "quadFit.hh"


namespace User
{
double* trajectoryTot(std::string theFile);
double* trajectoryToR(std::string theFile);
}

#endif
