/*
  This is a first pass at a program designed to perform a monte carlo
  simulation of the nuSol space probe's performance. The first 
  version is mean to be minimally functional, and will include wildly 
  inaccurate models.

  We generate a random number in [0,1]. The simulation will begin 
  stepping through the flight path until the probability of seeing a 
  neutrino event exceeds this random number. We will take this event 
  and give it energy and direction. Then we will smear these values 
  to represent our uncertainties. We apply an effeciency module and 
  then write this event to a ROOT tree. This continues until the 
  detector reaches the end of the flight path.
 */

//#include<random>
//#include<math.h>
#include"TROOT.h"
#include"neutrino.h"
/*
Double_t randomThreshold()
{
   double lower_bound = 0;
   double upper_bound = 1;
   std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
   std::default_random_engine re;
   double a_random_double = unif(re);

   return a_random_double;
   }
*/





namespace kUnits
{
  // Length (area, volume) Units.
  const Double_t meter = 1;
  const Double_t m = meter;

  const Double_t centimeter = meter/100;
  const Double_t cm = meter/100;

  const Double_t millimeter = meter/1000;
  const Double_t mm = meter/1000;

  const Double_t kilometer = 1000*meter;
  const Double_t km = 1000*meter;

  const Double_t solarRadii = 695700*km;

  const Double_t AU = 149597870700*m;



  // Time Units
  const Double_t second = 1;
  const Double_t s = 1;

  const Double_t milisecond = s/1000;
  const Double_t ms = s/1000;

  const Double_t microsecond = ms/1000;
  const Double_t us = ms/1000;

  const Double_t nanosecond = us/1000;
  const Double_t ns = nanosecond;

  const Double_t minute = 60*s;
  const Double_t min = 60*s;

  const Double_t hour = 60*minute;
  const Double_t hr = 60*minute;

  const Double_t day = 24*hr;
  const Double_t days = 24*hr;
  const Double_t d = 24*hr;

  const Double_t week = 7*day;
  const Double_t weeks = week;
  const Double_t wk = week;

  const Double_t year = 365.24219*day;
  const Double_t yr = 365.24219*day;



  // Mass Units




  // Energy Units




  // Angle Units




}



// This subprogram steps through the probe's flight path in
// finite time intervals. At each interval the program takes 
// the old distance and inclination value, and changes it to
// new values. The time step is currently hard-coded at a 
// constant value, and should be changed to 

void flightPath(Double_t &time, Double_t &distanceFromSunCenter, Double_t &inclination){
  inclination = 0;
  //std::cout << "distanceFromSunCenter :" << distanceFromSunCenter << endl;
  distanceFromSunCenter = 1*kUnits::AU-(time/kUnits::week)*kUnits::solarRadii; // spirals/falls inward
  //std::cout << "distanceFromSunCenter :" << distanceFromSunCenter << endl;
  time=time+1*kUnits::hour; // time steps in hours.
}




/*
Double_t *fightPath(Double_t stepNumber) {
  Double_t* array[3];
  for(Int_t i=0; i<3;i++){
    array[i]=stepNumber+i;
  }
  return array;
}
*/




// This program takes the detector's position (radius, inclination)
// and returns a double holding the cross section of a neutrino
// interaction. We want to code softly enough that the function
// can run in double-pulse and single pulse modes. 

Double_t neutrinoSignal (Double_t radius, Double_t inclination){
  return (1/TMath::Power( (radius/(1*kUnits::solarRadii)), 2));// Simple r^-2 law
}


// This program takes the detector's position (radius, inclination)
// and returns a double holding the cross section of a solar wind 
// background interaction. We want to code softly enough that the 
// function can run in double-pulse and single pulse modes. 
Double_t solarBackground (Double_t radius, Double_t inclination){
  return 1/TMath::Power( (radius/(1*kUnits::solarRadii)), 2);// Simple r^-2 law
}




// This program takes the detector's position (radius, inclination)
// and returns a double holding the cross section of a cosmic ray 
// background interaction. We want to code softly enough that the 
// function can run in double-pulse and single pulse modes. 
Double_t cosmicBackground (Double_t radius, Double_t inclination){
  return TMath::Exp(-1*kUnits::solarRadii/(radius)); // Goes to 0 at 0, and 1 far from the sun
}




// This program returns a probability of a radiological decay within a given timespan
Double_t radiologicalBackground (){
  return 123456789; // Goes to 0 at 0, and 1 far from the sun
}


// This is the main program
void nuSolPerformance(){
  // These are values that the main program manipulates
  Double_t time, radius, inclination=0;

  do{
    std::cout << "Before flightPath, time is " << (time/kUnits::days) << " d, radius is " << (radius/kUnits::solarRadii) << " RSol, and inclination is " << inclination <<endl;

    flightPath(time,radius,inclination);

    // std::cout << "After flightPath, time is " << time << " radius is " << radius << ", and inclination is " << inclination <<endl;

    
    bool cosmicBack = true;
    if(cosmicBack){  
      Double_t a =cosmicBackground(radius,1);
      std::cout << "The Cosmic Background of " << (radius/kUnits::solarRadii)  << " solarRadii from the center of the sun is : " << a << endl;
    }
  
    bool solarBack = true;
    if(solarBack){  
      Double_t b=solarBackground(radius,1);
      std::cout << "The Solar Background of " << (radius/kUnits::solarRadii)  << " solarRadii from the center of the sun is : " << b << endl;
    }
  
    bool neutrinos = true;
    if(neutrinos){  
      Double_t c = neutrinoSignal(radius,1);
      std::cout << "The Neutrino Signal of " << (radius/kUnits::solarRadii)  << " solarRadii from the center of the sun is : " << c << endl;
    }


    std::cout<<endl;
  }
  while(time<5*kUnits::year);

}


