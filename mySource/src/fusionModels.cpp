#include "fusionModels.hh"

namespace User
{
  void shellFusionModel(double radius, User::neutrinoData *theData){
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dist(0, 2*M_PI);
    theData -> neutrinoStartR = radius;
    theData -> neutrinoStartPhi = dist(gen);
    theData -> neutrinoStartTheta = dist(gen)/2;
  }
  
  
}
