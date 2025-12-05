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


  Version history (Should ahve been doing this more)
  - Jan 19: Exits
  Jan 19 2021 : Removed /steradian from the neutrino rates. (units in graph they were taken from did not have them. Oops)
  Feb 2       : Removing code into sub-files: kUnits,

  
 */

//#include<random>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>

#include <TROOT.h>
#include <TFile.h>
#include <TH1D.h>
#include <TRandom.h>
#include "kUnits.hh"
#include "quadFit.hh"
using namespace std;





// FLIGHT PATH STARTS -----------------------------------------


// This set of programs takes an input file of solar radii by
// day and outputs a quadratically interpolated value when
// handed a desired time of flight.
// -----------------------------------------------------
// DEPRICATED!!!!
// This subprogram steps through the probe's flight path in
// finite time intervals. At each interval the program takes 
// the old distance and inclination value, and changes it to
// new values. The time step is currently hard-coded at a 
// constant value, and should be changed to 
// ----------------------------------------------------





/*
// get number of lines  
int number_of_lines (std::string theFile) {
  // char* intermediate = theFile;
  std::ifstream intermediate (theFile);
  //data.open (intermediate);
  cout << "I have opened the file \"" << theFile << "\"\n\n";
  std::string line;
  int nLines=0;
  while (std::getline(intermediate, line))
    ++nLines;
  std::cout << "Number of lines in text file: " << nLines << endl;
  return nLines;
}
  



// This subprogram interpolates parabolas from the file we hand it
// We choose parabolas to account for the changing sign when our 
// orbit hits an edge. This function takes off two days.


double *johnTheInterpolator(std::string theFile){
  ifstream data;
  // get number of lines  
  int nLines = number_of_lines(theFile);

  double datum;
  data.open (theFile);
  cout << "I have opened the file \"" << theFile << " for the second time\"\n\n";
  int count=0;
  double y[nLines];// intermediate values  y(0 days), y(1 day), ...
  while(data >> datum){
    cout << "Line " << count << " gives datum " << datum <<endl;
    y[count]=datum;
    count++;
  }
  cout << "I've collected the data set.\n";
  for(int i = 0; i<nLines; i++){
    //cout << intermediate[i] << endl;
  }

  // calculate a, b, c for at^2+bt+c in order to fit the radius as
  // a function of time for day 
  double *theValues = new  double[3*(nLines-2)];// return the values as one long array
  for(int i=1; i<(nLines-1);i++){
    // values found from calculating (a,b,c)
    // =inverse( (1,x(i-1),x^2(i-1); 1,x(i),x^2(i)); 1,x(i+1),x^2(i+1) )
    // (inner product) (c,b,a)^T
    
    double a = y[i-1]/2-y[i]+y[i+1]/2;
    cout << "a[" << i << "] is " << a << endl;
      // ( 2*y[i]-y[i-1]-y[i+1] ) // numerator
      // / ( -2 ); // denominator
    double b = ( (-2*i-1)*y[i-1] + 4*i*y[i] + (-2*i+1)*y[i+1] )/2;
      //( 2*y[i-1] - y[i] - y[i+1] ) / ( -3*i ) 
      //	- a * ( -6*i-1 ) / ( -3 );
    double c = ( (i*i+i)*y[i-1] + (2-2*i*i)*y[i] + (i*i-i)*y[i+1] )/2;
      //( y[i-1] + y[i] + y[i+1] + a*( 3*i+2 ) + b*( 3*i ) ) / 3;
    cout << "a[" << i << "] is " << a << ". b[" << i << "] is " << b << ". c[" << i << "] is " << c  << endl;
    
    // data for returning
    theValues [3*(i-1)]=c; // x^0 constant
    theValues [3*(i-1)+1]=b; // x^1 constant
    theValues [3*(i-1)+2]=a; // x^2 constant
  }

  return theValues;
}



double* elliptical (double closest, double furthest){
  cout << "I am " << furthest/kUnits::AU << " AU at my furthest.\n";
  cout << "I am " << closest/kUnits::solarRadii << " Rsol at my closest.\n";
  
  cout << "Furthest - closext = " << ((furthest-closest)/kUnits::solarRadii) << endl;
  cout << "Furthest + closext = " << ((furthest+closest)/kUnits::solarRadii) << endl;
  double e = (furthest-closest)/(furthest+closest);
  cout << "The eecentricity is " << e << ".\n";
  double p = furthest*(1-e);
  double a = (closest+furthest)/2;// semi-major axis
  double H = furthest*sqrt(kUnits::G*kUnits::mSun* ( 2/furthest - 1/a ) ) ;//  specific angular momentum
  const double timeStep = M_PI/10*kUnits::sec;  cout << "The specific angular momentum is " << H << endl;

  double angle = -M_PI;// angle at the start; at furthest approach
  double radius =  p/(1 + e*cos(angle) );// r(theta)
  double timeMax = 2*M_PI*sqrt( pow(a,3) / ( kUnits::G*kUnits::mSun ) );
  cout << "The orbital period is " << timeMax/kUnits::day << " days.\n";
  double nSteps = timeMax/timeStep;
  int intSteps = ceil(nSteps);
  long long int cells = ceil(timeMax/kUnits::min);
  //double theValues[cells];
  double theValues[cells];
  int sanityCount = 0;


  cout << "I'm about to find the position/radii/etc for the ellipse uising " << intSteps << " steps. I am going to put them into " << cells << " cells.\n\n";
  for (long long int i = 0; i < intSteps; i++){
    int cellNum = ceil(i*timeStep/kUnits::min);    
    //    cout << "day = " << i*timeStep/kUnits::days << endl;
    if (i%10000==0){
    cout << "r = " << radius/kUnits::solarRadii << endl;
    theValues[cellNum]=radius/kUnits::solarRadii;// in units of solar radii
    //cout << "I have written " << theValues[cellNum] << " to the array.\n";
    }
    double dAngle = H/(radius)/(radius);// dtheta/dt
    angle = angle + dAngle*timeStep;// update angle for dt
    radius = p/(1 + e*cos(angle) ); // update r for new angle



    
    //cout << "r = " << radius/kUnits::solarRadii << endl;
    //cout << "I'm filling cell number " << cellNum << " with " << radius/kUnits::solarRadii << endl;
    theValues[cellNum]=radius/kUnits::solarRadii;// in units of solar radii
    
    // Data updating
    
  }
  cout << "There are " << sizeof(theValues)/sizeof(theValues[0]) << " or " << cells << " values in the radius array.\n";

  // Check the radii
  for (int i=0;i<cells;i++){
    if(i%50000==0){
      cout << "The array value for the " << i << "th radius is r = " << theValues[i] << endl;
    }
  }


  cout << "\n\nI'm about to find the fit parameters for the ellipse, which should hold " << (3*(cells-1))  << "values.\n\n";
  double* theReturns = new double[(3*(cells-1))];// return the values as one long array
  // this type of declaration shuld fix a memory issue. We keep theReturns around until deleted


  for(long long int i=1; i<(cells-1);i++){
    // values found from calculating (a,b,c)
    // =inverse( (1,x(i-1),x^2(i-1); 1,x(i),x^2(i)); 1,x(i+1),x^2(i+1) )
    // (inner product) (c,b,a)^T
    double a = theValues[i-1]/2-theValues[i]+theValues[i+1]/2;
    // ( 2*theValues[i]-theValues[i-1]-theValues[i+1] ) // numerator
    // / ( -2 ); // denominator
    double b = ( (-2*i-1)*theValues[i-1] + 4*i*theValues[i] + (-2*i+1)*theValues[i+1] )/2;
    //( 2*theValues[i-1] - theValues[i] - theValues[i+1] ) / ( -3*i ) 
    //      - a * ( -6*i-1 ) / ( -3 );
    double c = ( (i*i+i)*theValues[i-1] + (2-2*i*i)*theValues[i] + (i*i-i)*theValues[i+1] )/2;
    
    if( i%100==0){
    cout << "c = " << c <<".\n";
    }

    if (c>1000000){
      cout << "c is more than a million... somehow. The iteration number is "
	   << i << ", and the cell number is " << (3*(i-1)) << ".\nFor good measure, b = " 
	   << b <<", and a = " << a << endl;
	break;
    }  
    if (c<-1000000){
      cout << "c is less than negative 1 million... somehow. The cell number is "	   << i << ".\n";
	break;
    }  
    //( theValues[i-1] + theValues[i] + theValues[i+1] + a*( 3*i+2 ) + b*( 3*i ) ) / 3;
    //cout << "a[" << i-1 << "] is " << a << ". b[" << i-1 << "] is " << b << ". c[" << i-1 << "] is " << c  << endl;
    //cout << "Hence the radius at day " << i << " is " << (a*i*i + b*i + c    ) << ".\n\n";

    // data for returning
    theReturns [3*(i-1)]=c; // x^0 constant
    theReturns [3*(i-1)+1]=b; // x^1 constant
    theReturns [3*(i-1)+2]=a; // x^2 constant
  }
  cout << "c[0], b[0], and a[0], are " << theReturns[0] << ", " << theReturns[1] << ", " << theReturns[2] << " respectively.\n\n";

  cout << "elliptical is done iterating.\n";
  
  return theReturns;
}




double positionDays(double time, double* quadConst, double maxTime){
  //cout << "Position function has been called.\n";
  double currentRadius=0;
  int day = floor(time/kUnits::day);
  // cout is slow
    //cout << "The time is " << (time/kUnits::hr) << "hours in, and we are " 
      // << day << " days in.\n";
  //
  if(time < (1*kUnits::day)){// first day from later data
    currentRadius = *(quadConst)+time*(*( quadConst+1 ) )/kUnits::day 
      + time * time * (*( quadConst+2 ) ) /kUnits::day /kUnits::day;
  }
  else if (time > (maxTime-kUnits::day)){ // last day from earlier data
    currentRadius = *(quadConst + 3*( day-2 ) )
      + *( quadConst+1  + 3*( day-2 ) ) * time /kUnits::day 
      + *( quadConst+2  + 3*( day-2 ) ) * time * time /kUnits::day /kUnits::day;    
  }
  else{
    currentRadius = *(quadConst + 3*( day-1 ) )
      + *( quadConst+1  + 3*( day-1 ) ) * time /kUnits::day 
      + *( quadConst+2  + 3*( day-1 ) ) * time * time /kUnits::day /kUnits::day;
  }
  //cout << "The current radius is " << currentRadius << endl;
  return currentRadius*kUnits::solarRadii; // return the radius in solar radii
}

// hours version
double positionHours(double time, double* quadConst, double maxTime){
  //cout << "Position function has been called.\n";
  double currentRadius=0;
  int day = floor(time/kUnits::day);
  //cout << "The time is " << (time/kUnits::hr) << "hours in, and we are " 
  //     << day << " days in.\n";
  if(time < (1*kUnits::hr)){
    currentRadius = *(quadConst)+time*(*( quadConst+1 ) )/kUnits::hr 
      + time * time * (*( quadConst+2 ) ) /kUnits::hr /kUnits::hr;
  }
  else if (time > (maxTime - kUnits::hr)){ // last day from earlier data
    currentRadius = *(quadConst + 3*( day-2 ) )
      + *( quadConst+1  + 3*( day-2 ) ) * time /kUnits::hr 
      + *( quadConst+2  + 3*( day-2 ) ) * time * time /kUnits::hr /kUnits::hr;    
    cout << "This is the last time set. a = " 
	 << *( quadConst+2  + 3*( day-2 ) ) << ". b = " 
	 << *( quadConst+1  + 3*( day-2 ) ) << ". c = " 
	 << *( quadConst+2  + 3*( day-2 ) ) << ". t = " << time 
	 << ". Hence the radius is " 
	 << (*(quadConst + 3*( day-2 ) )
	     + *( quadConst+1  + 3*( day-2 ) ) * time /kUnits::hr 
	     + *( quadConst+2  + 3*( day-2 ) ) * time * time /kUnits::hr /kUnits::hr)
	 << ".\n";
  }
  else{
    int day = floor(time/kUnits::hr);
    // couts are slow
    //cout << "a = " << *( quadConst+2  + 3*( day-1 ) ) << ". b = "
//	 << *( quadConst+1  + 3*( day-1 ) ) << ". c = "
//	 << *(quadConst + 3*( day-1 ) ) << ".\n";
    

    currentRadius = *(quadConst + 3*( day-1 ) )
      + *( quadConst+1  + 3*( day-1 ) ) * time /kUnits::hr 
      + *( quadConst+2  + 3*( day-1 ) ) * time * time /kUnits::hr /kUnits::hr;
    
  }
  //cout << "The current radius is " << currentRadius << endl;
  return currentRadius*kUnits::solarRadii; // return the radius in solar radii
}

// minutes version
double positionMinutes(double time, double* quadConst, double maxTime){
  //cout << "Position function has been called.\n";
  double currentRadius=0;
  int day = floor(time/kUnits::day);
  //cout << "The time is " << (time/kUnits::hr) << "hours in, and we are " 
  //     << day << " days in.\n";
  if(time < (1*kUnits::min)){
    currentRadius = *(quadConst)+time*(*( quadConst+1 ) )/kUnits::min 
      + time * time * (*( quadConst+2 ) ) /kUnits::min /kUnits::min;
  }
  else if (time > (maxTime - kUnits::hr)){ // last day from earlier data
    currentRadius = *(quadConst + 3*( day-2 ) )
      + *( quadConst+1  + 3*( day-2 ) ) * time /kUnits::minute 
      + *( quadConst+2  + 3*( day-2 ) ) * time * time /kUnits::minute /kUnits::min;    
    
    cout << "This is the last time set. a = " << *( quadConst+2  + 3*( day-2 ) ) 
	 << ". b = " << *( quadConst+1  + 3*( day-2 ) ) << ". c = " 
	 << *( quadConst+2  + 3*( day-2 ) ) << ". t = " << time 
	 << ". Hence the radius is " << (*(quadConst + 3*( day-2 ) )
					 + *( quadConst+1  + 3*( day-2 ) ) * time /kUnits::minute 
					 + *( quadConst+2  + 3*( day-2 ) ) * time * time /kUnits::minute /kUnits::min)
	 << ".\n";
  }
  else{
    int day = floor(time/kUnits::min);
    // couts are slow
      //cout << "a = " << *( quadConst+2  + 3*( day-1 ) ) << ". b = "
        //   << *( quadConst+1  + 3*( day-1 ) ) << ". c = "
	  // << *(quadConst + 3*( day-1 ) ) << ".\n";
    //
    currentRadius = *(quadConst + 3*( day-1 ) )
      + *( quadConst+1  + 3*( day-1 ) ) * time /kUnits::min 
      + *( quadConst+2  + 3*( day-1 ) ) * time * time /kUnits::min /kUnits::min;
    
  }
  // cout << "The current radius is " << currentRadius << endl;
  return currentRadius*kUnits::solarRadii; // return the radius in solar radii
}
*/

/*
void flightPath(double &time, double &distanceFromSunCenter, double &inclination){
  inclination = 0;
  //std::cout << "distanceFromSunCenter :" << distanceFromSunCenter << endl;
  //distanceFromSunCenter = 1*kUnits::AU-(time/kUnits::week)*kUnits::solarRadii; // spirals/falls inward
  //oscillates from 7 to 147 solar radii with period 2*pi hours
  // distanceFromSunCenter =14*( 1+0.5*sin( time/ (kUnits::hour) ) )*kUnits::solarRadii;
  distanceFromSunCenter = 3*kUnits::solarRadii;

  //std::cout << "distanceFromSunCenter :" << distanceFromSunCenter << endl;
  time=time+1*timeStep; // time steps in hours.
}

*/

// FLIGHT PATH ENDS ------------------------------------------





// NEUTRINOS START ---------------------------------------------



// helper function for radius to energy in focus calcs
double rToE(double radius){
  double startRad = 19*kUnits::AU;
  double endRad = 25*kUnits::AU;
  double startEner = 0.1*kUnits::MeV;
  double endEner = 20*kUnits::MeV;
  double slope = (endEner-startEner)/(endRad-startRad);
  double intercept = endEner - slope*endRad;
  double energy = slope*radius - intercept;
  return energy;
}

// NEUTRINOS BY PROCESS AND ENERGY


double hepFlux (double eMin, double eMax){
  eMin = eMin/kUnits::MeV;// fix units for integrating
  eMax = eMax/kUnits::MeV;
  
  // 5 regions with 5 fits
  double gMin = 0.2332; // Gallium interaction threshold determines region 1 minimum
  double p1ConstPoint = 7.6; // Switch from power1 to constant
  double constP2Point = 10.5; // Switch from constant to power2
  double p2P3Point = 14; // Switch from power2 to power3
  double p3ZeroPoint = 15.5; // switch from power3 to 0

  // values from IntegratingNeutrinoProbability.ods
  double const1 = 26.31;
  double pow1 = 1.618;
  double constant = 700;
  double const2 = 712465;
  double pow2 = -2.945;
  double const3 = 2.295e17;
  double pow3 = -12.99;

  // Values for integral evaluation
  // if eMin>max => max
  // if max>eMin>min => eMin
  // if min>eMax => min
  double pow1Min = min( p1ConstPoint , max( eMin , gMin ) );
  double constMin = min( constP2Point , max( eMin , p1ConstPoint ) );
  double pow2Min = min( p2P3Point , max( eMin , constP2Point ) );
  double pow3Min = min( p3ZeroPoint , max( eMin , p2P3Point ) );
  
  // idential except eMin->eMax
  double pow1Max = min( p1ConstPoint , max( eMax , gMin ) );
  double constMax = min( constP2Point , max( eMax , p1ConstPoint ) );
  double pow2Max = min( p2P3Point , max( eMax , constP2Point ) );
  double pow3Max = min( p3ZeroPoint , max( eMax , p2P3Point ) );

  double integral = 0;// ininitalize 
  integral += const1*( pow(pow1Max,pow1+1) - pow(pow1Min,pow1+1) ) / (pow1+1); // Power integral
  integral += const2*( pow(pow2Max,pow2+1) - pow(pow2Min,pow2+1) ) / (pow2+1);
  integral += const3*( pow(pow3Max,pow3+1) - pow(pow3Min,pow3+1) ) / (pow3+1);
  integral += constant*(constMax-constMin);// constant integral
  
  return integral/(kUnits::cm)/(kUnits::cm)/(kUnits::s);// per second per cm^2
}

double B8Flux (double eMin, double eMax){
  eMin = eMin/kUnits::MeV;// fix units for integrating
  eMax = eMax/kUnits::MeV;
  
  // 5 regions with 5 fits
  double gMin = 0.2332; // Gallium interaction threshold determines region 1 minimum
  double p1ConstPoint = 5.05; // Switch from power1 to constant
  double constP2Point = 8; // Switch from constant to power2
  double p2P3Point = 10.05; // Switch from power2 to power3
  double p3ZeroPoint = 12; // switch from power3 to 0

  // values from IntegratingNeutrinoProbability.ods
  double const1 = 408558;
  double pow1 = 0.3325;
  double constant = 7e5;
  double const2 = 2455031377;
  double pow2 = -3.863;
  double const3 = 2.492e15;
  double pow3 = -9.913;

  // Values for integral evaluation
  // if eMin>max => max
  // if max>eMin>min => eMin
  // if min>eMax => min
  double pow1Min = min( p1ConstPoint , max( eMin , gMin ) );
  double constMin = min( constP2Point , max( eMin , p1ConstPoint ) );
  double pow2Min = min( p2P3Point , max( eMin , constP2Point ) );
  double pow3Min = min( p3ZeroPoint , max( eMin , p2P3Point ) );
  
  // idential except eMin->eMax
  double pow1Max = min( p1ConstPoint , max( eMax , gMin ) );
  double constMax = min( constP2Point , max( eMax , p1ConstPoint ) );
  double pow2Max = min( p2P3Point , max( eMax , constP2Point ) );
  double pow3Max = min( p3ZeroPoint , max( eMax , p2P3Point ) );

  double integral = 0;// ininitalize 
  integral += const1*( pow(pow1Max,pow1+1) - pow(pow1Min,pow1+1) ) / (pow1+1); // Power integral
  integral += const2*( pow(pow2Max,pow2+1) - pow(pow2Min,pow2+1) ) / (pow2+1);
  integral += const3*( pow(pow3Max,pow3+1) - pow(pow3Min,pow3+1) ) / (pow3+1);
  integral += constant*(constMax-constMin);// constant integral
  
  return integral/(kUnits::cm)/(kUnits::cm)/(kUnits::s);// per second per cm^2
}

double F17Flux (double eMin, double eMax){
  eMin = eMin/kUnits::MeV;// fix units for integrating
  eMax = eMax/kUnits::MeV;
  
  // 5 regions with 5 fits
  double gMin = 0.2332; // Gallium interaction threshold determines region 1 minimum
  double p1ConstPoint = 0.95; // Switch from power1 to constant
  double constP2Point = 1.06; // Switch from constant to power2
  double p2P3Point = 1.4; // Switch from power2 to power3
  double p3ZeroPoint = 1.7; // switch from power3 to 0

  // values from IntegratingNeutrinoProbability.ods
  double const1 = 6499028;
  double pow1 = 1.558;
  double constant = 6e6;
  double const2 = 6937470;
  double pow2 = -2.4915;
  double const3 = 98537561;
  double pow3 = -10.38;

  // Values for integral evaluation
  // if eMin>max => max
  // if max>eMin>min => eMin
  // if min>eMax => min
  double pow1Min = min( p1ConstPoint , max( eMin , gMin ) );
  double constMin = min( constP2Point , max( eMin , p1ConstPoint ) );
  double pow2Min = min( p2P3Point , max( eMin , constP2Point ) );
  double pow3Min = min( p3ZeroPoint , max( eMin , p2P3Point ) );
  
  // idential except eMin->eMax
  double pow1Max = min( p1ConstPoint , max( eMax , gMin ) );
  double constMax = min( constP2Point , max( eMax , p1ConstPoint ) );
  double pow2Max = min( p2P3Point , max( eMax , constP2Point ) );
  double pow3Max = min( p3ZeroPoint , max( eMax , p2P3Point ) );

  double integral = 0;// ininitalize 
  integral += const1*( pow(pow1Max,pow1+1) - pow(pow1Min,pow1+1) ) / (pow1+1); // Power integral
  integral += const2*( pow(pow2Max,pow2+1) - pow(pow2Min,pow2+1) ) / (pow2+1);
  integral += const3*( pow(pow3Max,pow3+1) - pow(pow3Min,pow3+1) ) / (pow3+1);
  integral += constant*(constMax-constMin);// constant integral
  
  return integral/(kUnits::cm)/(kUnits::cm)/(kUnits::s);// per second per cm^2
}

double O15Flux (double eMin, double eMax){
  eMin = eMin/kUnits::MeV;// fix units for integrating
  eMax = eMax/kUnits::MeV;
  
  // 5 regions with 5 fits
  double gMin = 0.2332; // Gallium interaction threshold determines region 1 minimum
  double p1ConstPoint = 0.94; // Switch from power1 to constant
  double constP2Point = 1.05; // Switch from constant to power2
  double p2P3Point = 1.3; // Switch from power2 to power3
  double p3ZeroPoint = 1.6; // switch from power3 to 0

  // values from IntegratingNeutrinoProbability.ods
  double const1 = 219398872;
  double pow1 = 1.496;
  double constant = 2e8;
  double const2 = 229267829;
  double pow2 = -2.799;
  double const3 = 948164019;
  double pow3 = -8.201;

  // Values for integral evaluation
  // if eMin>max => max
  // if max>eMin>min => eMin
  // if min>eMax => min
  double pow1Min = min( p1ConstPoint , max( eMin , gMin ) );
  double constMin = min( constP2Point , max( eMin , p1ConstPoint ) );
  double pow2Min = min( p2P3Point , max( eMin , constP2Point ) );
  double pow3Min = min( p3ZeroPoint , max( eMin , p2P3Point ) );
  
  // idential except eMin->eMax
  double pow1Max = min( p1ConstPoint , max( eMax , gMin ) );
  double constMax = min( constP2Point , max( eMax , p1ConstPoint ) );
  double pow2Max = min( p2P3Point , max( eMax , constP2Point ) );
  double pow3Max = min( p3ZeroPoint , max( eMax , p2P3Point ) );

  double integral = 0;// ininitalize 
  integral += const1*( pow(pow1Max,pow1+1) - pow(pow1Min,pow1+1) ) / (pow1+1); // Power integral
  integral += const2*( pow(pow2Max,pow2+1) - pow(pow2Min,pow2+1) ) / (pow2+1);
  integral += const3*( pow(pow3Max,pow3+1) - pow(pow3Min,pow3+1) ) / (pow3+1);
  integral += constant*(constMax-constMin);// constant integral
  
  return integral/(kUnits::cm)/(kUnits::cm)/(kUnits::s);// per second per cm^2
}

double N13Flux (double eMin, double eMax){
  eMin = eMin/kUnits::MeV;// fix units for integrating
  eMax = eMax/kUnits::MeV;
  
  // 5 regions with 5 fits
  double gMin = 0.2332; // Gallium interaction threshold determines region 1 minimum
  double p1ConstPoint = 0.7; // Switch from power1 to constant
  double constP2Point = 0.84; // Switch from constant to power2
  double p2P3Point = 1; // Switch from power2 to power3
  double p3ZeroPoint = 1.1; // switch from power3 to 0

  // values from IntegratingNeutrinoProbability.ods
  double const1 = 695177565;
  double pow1 = 1.48;
  double constant = 4.1e8;
  double const2 = 3.4e8;
  double pow2 = -1.074;
  double const3 = 3.4e8;
  double pow3 = -13.95;

  // Values for integral evaluation
  // if eMin>max => max
  // if max>eMin>min => eMin
  // if min>eMax => min
  double pow1Min = min( p1ConstPoint , max( eMin , gMin ) );
  double constMin = min( constP2Point , max( eMin , p1ConstPoint ) );
  double pow2Min = min( p2P3Point , max( eMin , constP2Point ) );
  double pow3Min = min( p3ZeroPoint , max( eMin , p2P3Point ) );
  
  // idential except eMin->eMax
  double pow1Max = min( p1ConstPoint , max( eMax , gMin ) );
  double constMax = min( constP2Point , max( eMax , p1ConstPoint ) );
  double pow2Max = min( p2P3Point , max( eMax , constP2Point ) );
  double pow3Max = min( p3ZeroPoint , max( eMax , p2P3Point ) );

  double integral = 0;// ininitalize 
  integral += const1*( pow(pow1Max,pow1+1) - pow(pow1Min,pow1+1) ) / (pow1+1); // Power integral
  integral += const2*( pow(pow2Max,pow2+1) - pow(pow2Min,pow2+1) ) / (pow2+1);
  integral += const3*( pow(pow3Max,pow3+1) - pow(pow3Min,pow3+1) ) / (pow3+1);
  integral += constant*(constMax-constMin);// constant integral
  
  return integral/(kUnits::cm)/(kUnits::cm)/(kUnits::s);// per second per cm^2
}

double ppFlux (double eMin, double eMax){
  eMin = eMin/kUnits::MeV;// fix units for integrating
  eMax = eMax/kUnits::MeV;
  
  // 5 regions with 5 fits
  double gMin = 0.2332; // Gallium interaction threshold determines region 1 minimum
  double p1ConstPoint = 0.25; // Switch from power1 to constant
  double constP2Point = 0.32; // Switch from constant to power2
  double p2P3Point = 0.39; // Switch from power2 to power3
  double p3ZeroPoint = 0.405; // switch from power3 to 0

  // values from IntegratingNeutrinoProbability.ods
  double const1 = 1733245360905;
  double pow1 = 1.54;
  double constant = 2.05e11;
  double const2 = 96924426898;
  double pow2 = -0.6574;
  double const3 = 0.2245;
  double pow3 = -29.11;

  // Values for integral evaluation
  // if eMin>max => max
  // if max>eMin>min => eMin
  // if min>eMax => min
  double pow1Min = min( p1ConstPoint , max( eMin , gMin ) );
  double constMin = min( constP2Point , max( eMin , p1ConstPoint ) );
  double pow2Min = min( p2P3Point , max( eMin , constP2Point ) );
  double pow3Min = min( p3ZeroPoint , max( eMin , p2P3Point ) );
  
  // idential except eMin->eMax
  double pow1Max = min( p1ConstPoint , max( eMax , gMin ) );
  double constMax = min( constP2Point , max( eMax , p1ConstPoint ) );
  double pow2Max = min( p2P3Point , max( eMax , constP2Point ) );
  double pow3Max = min( p3ZeroPoint , max( eMax , p2P3Point ) );

  double integral = 0;// ininitalize 
  integral += const1*( pow(pow1Max,pow1+1) - pow(pow1Min,pow1+1) ) / (pow1+1); // Power integral
  integral += const2*( pow(pow2Max,pow2+1) - pow(pow2Min,pow2+1) ) / (pow2+1);
  integral += const3*( pow(pow3Max,pow3+1) - pow(pow3Min,pow3+1) ) / (pow3+1);
  integral += constant*(constMax-constMin);// constant integral
  
  //cout << "The pp integral flux for " << eMin << " to " << eMax << " in MeV is " << integral << endl;
  
  return integral/(kUnits::cm)/(kUnits::cm)/(kUnits::s);// per second per cm^2
}

double pepFlux(double eMin, double eMax){
  double eDelta = 1.3*kUnits::MeV; //placeholder
  if ( (eMin < eDelta) && (eMax > eDelta)){
    return 1.3e8/(kUnits::cm)/(kUnits::cm)/(kUnits::s);// per second per cm^2
  }
  else{
    return 0;
  } 
}

double Be7Flux(double eMin, double eMax){
  double eDelta1 = 0.38*kUnits::MeV; //placeholder
  double eDelta2 = 0.87*kUnits::MeV; //placeholder
  if ( (eMin < eDelta1) && (eMax > eDelta1)){
    return 5e8/(kUnits::cm)/(kUnits::cm)/(kUnits::s);// per second per cm^2
  }
  else if ( (eMin < eDelta2) && (eMax > eDelta2)){
    return 4.1e9/(kUnits::cm)/(kUnits::cm)/(kUnits::s);// per second per cm^2
  }
  else{
    return 0;
  } 
}

// Cross sections for processes

double ppCrossSec(double eMin, double eMax){
  return 11.72e-46*kUnits::cm*kUnits::cm;
}

double pepCrossSec(double eMin, double eMax){
  return 204e-46*kUnits::cm*kUnits::cm;
}

double Be7CrossSec(double eMin, double eMax){
  return 71.7e-46*kUnits::cm*kUnits::cm;
}

double N13CrossSec(double eMin, double eMax){
  return 60.4e-46*kUnits::cm*kUnits::cm;
}

double O15CrossSec(double eMin, double eMax){
  return 113.7e-46*kUnits::cm*kUnits::cm;
}

double F17CrossSec(double eMin, double eMax){
  return 113.9e-46*kUnits::cm*kUnits::cm;
}


double B8CrossSec(double eMin, double eMax){
  return 2.4e-42*kUnits::cm*kUnits::cm;
}

double hepCrossSec(double eMin, double eMax){
  return 7.14e-42*kUnits::cm*kUnits::cm;
}

double kgGallium = 200;
double galliumAtomNumber = kgGallium * 8.637e24;//

// This program takes the detector's position (radius, inclination)
// and returns a double holding the flux of  neutrinos at that
// point. We want to code softly enough that the function
// can run in double-pulse and single pulse modes. 

double neutrinoSignal (double radius, double inclination, double eMin, double eMax){
  
  // number of interactions per second per atom at 1 AU. 
  // Calculated from effective cross sections and solar
  // neutrino flux by source graph. 
  double neutrinoInteractionRate = 0;//1.05957e-34;

  // I forgot to add a for loop to integrate over the energy range. I want to do it logrythmically

  // Integral rate
  neutrinoInteractionRate += ppCrossSec(eMin, eMax) * ppFlux(eMin, eMax);
  neutrinoInteractionRate += pepCrossSec(eMin, eMax) * pepFlux(eMin, eMax);
  neutrinoInteractionRate += Be7CrossSec(eMin, eMax) * Be7Flux(eMin, eMax);
  neutrinoInteractionRate += N13CrossSec(eMin, eMax) * N13Flux(eMin, eMax);
  neutrinoInteractionRate += O15CrossSec(eMin, eMax) * O15Flux(eMin, eMax);
  neutrinoInteractionRate += F17CrossSec(eMin, eMax) * F17Flux(eMin, eMax);
  neutrinoInteractionRate += B8CrossSec(eMin, eMax) * B8Flux(eMin, eMax);
  neutrinoInteractionRate += hepCrossSec(eMin, eMax) * hepFlux(eMin, eMax);
  /* 
 cout << "The pp rate is  " << ppCrossSec(eMin, eMax) * ppFlux(eMin, eMax) << " and flux is " << ppFlux(eMin, eMax)  <<endl;
  cout << "The pep rate is " << pepCrossSec(eMin, eMax) * pepFlux(eMin, eMax) << " and flux is " << pepFlux(eMin, eMax)  << endl;
  cout << "The Be7 rate is " << Be7CrossSec(eMin, eMax) * Be7Flux(eMin, eMax) << " and flux is " << Be7Flux(eMin, eMax)  << endl;
  cout << "The N13 rate is " << N13CrossSec(eMin, eMax) * N13Flux(eMin, eMax) << " and flux is " << N13Flux(eMin, eMax)  << endl;
  cout << "The O15 rate is " << O15CrossSec(eMin, eMax) * O15Flux(eMin, eMax) << " and flux is " << O15Flux(eMin, eMax)  << endl;
  cout << "The F17 rate is " << F17CrossSec(eMin, eMax) * F17Flux(eMin, eMax) << " and flux is " << F17Flux(eMin, eMax)  << endl;
  cout << "The B8 rate is  " << B8CrossSec(eMin, eMax) * B8Flux(eMin, eMax) << " and flux is " << B8Flux(eMin, eMax)  << endl;
  cout << "The hep rate is " << hepCrossSec(eMin, eMax) * hepFlux(eMin, eMax) << " and flux is " << hepFlux(eMin, eMax) << endl;
  */

  // Fiat interaction rate
  //neutrinoInteractionRate = 1.05957*pow(10,-34);

  //cout << "The total interaction rate is " << neutrinoInteractionRate << "per gallium atom.\n";

  // find the rate per kg of gallium
  neutrinoInteractionRate = neutrinoInteractionRate*galliumAtomNumber;
  
  // cout << "The total interaction rate is " << neutrinoInteractionRate << "for the total gallium mass.\n";
  
  // find the rate at arbitrary radius
  neutrinoInteractionRate = neutrinoInteractionRate/( (radius/kUnits::AU)*(radius/kUnits::AU));
  //cout << "The total interaction rate is " << neutrinoInteractionRate << "at this particular radius.\n";
  
  return neutrinoInteractionRate;
}



// NEUTRINOS END ---------------------------------------------------------








// BACKGROUNDS START ----------------------------------------------------


// This function determines the neutron cross section

double neutronCrossSection (double radius, double inclination, double energy){
  return (10 *kUnits::barns);// approximately constant or smaller than this for all energies of interest.
}

// This program uses approximations from
// https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2011SW000732
// to find the neutron flux. Assuming each flare has approximately the 
// same overall flux, we can find the flux from one 



double neutronFlux (double radius, double inclination, double energy){
/*
  double intermediate;
  double mc2 = 939.565560*kUnits::MeV;
  // assume: mean lifetime =886 sec, mc^2=939.565560 MeV
  double survivalProb = exp( -(radius-1)/(886*kUnits::seconds*kUnits::c*sqrt( pow( ((energy+mc2)/mc2),2)-1 ) ) );


  if((energy/kUnits::MeV) < 1){
 
    intermediate =  (3*pow(10,-4)) / (kUnits::MeV);
  }  
  else if((energy/kUnits::MeV) < 70){
    intermediate =  ( pow(energy,-0.5830)*3*pow(10,-4) ) / (kUnits::MeV);
  }
  else{
    intermediate =  ( pow(energy,-1.862)*pow(10,-2.564) ) / (kUnits::MeV);
  }

  return intermediate*pow( (radius/kUnits::solarRadii) , -2) / (3*kUnits::hr)*pow(10,28);
*/
  return (1e-5 / pow(kUnits::cm,2))/kUnits::second/pow(radius/kUnits::solarRadii,2);// near zero outsie of solar flares?
}

// This program takes the detector's position (radius, inclination)
// and returns a double holding the cross section of a solar wind 
// background interaction. We want to code softly enough that the 
// function can run in double-pulse and single pulse modes.

double neutronTargetNumber = 1e27; // 10 times gallium target number for now
double protonTargetNumber = 1e28; // ~assuming the target is 400 pounds of tungsten

double solarBackground (double radius, double inclination, double detectorDiameter){
  // assuming the tungsten is 5 cm thick, the stopping power is 5*20
  // multiplied by the data at https://physics.nist.gov/cgi-bin/Star/ap_table.pl
  // This means we only need to worry about protons >~ 10 MeV
  // We'll assume the 10 MeV warning threshold from the GOES data,
  // 10 cm^-2 s^-1 sr^-1, as a basis for analysis.
  // If this is too high, we can look more in depth for solar maximum data.
  // https://www.swpc.noaa.gov/products/goes-proton-flux-dynamic-plot
  double area = M_PI*radius*radius;// r is about 0.3 meters for our detector
  // but there's also some small bitys of the sides that will be hit
  

  
  
  // Find the effective radius of the area completely covered by the detector
  double rUmbra = (kUnits::solarRadii-detectorDiameter/2) * 
    ( 1*kUnits::AU/radius + 1/(1-detectorDiameter/(2*kUnits::solarRadii)) - 1 );
  double aUmbra = M_PI*rUmbra*rUmbra;
  // And partially covered
  double rPenumbra = detectorDiameter*(1+2*kUnits::solarRadii/detectorDiameter)*
    (1*kUnits::AU/radius + (1/( 1 + 2*kUnits::solarRadii/detectorDiameter ) -1  ) );
  // assume scale factor linear from 0 at penumbra radius to 1 at umbra radius
  double aPenumbra = 2*M_PI*( (pow(rUmbra,3)- pow(rPenumbra,3)) / (3*(rPenumbra-rUmbra) ) 
			      + (1 + rUmbra / (rPenumbra-rUmbra) )
			      * ( (-pow(rUmbra,2) + pow(rPenumbra,2) ) / 2) );
  
  /*
    cout << "The current radius is " << radius/kUnits::solarRadii << " solar radii.\n";
    cout << "The radius of the umbra is " << rUmbra/kUnits::km << " km, and the area is " << aUmbra/(kUnits::km*kUnits::km) << " km^2.\n"; 
    cout << "The radius of the penumbra is " << rPenumbra/kUnits::km << " km, and the area is " << aPenumbra/(kUnits::km*kUnits::km) << " km^2.\n";
  */
  double angularView = M_PI;// quarter of view, probably alwasys smaller than this, but dependent on radius
                              // this might take the 1/r^2 depndence term in later.
  
  double earthValue = 10/kUnits::cm/kUnits::cm/kUnits::sec;// sr^-1
  double flux = earthValue;
  // for shadow method
  //flux = flux*(aUmbra+aPenumbra);
  
  // for 1/r^2 method
  earthValue *(kUnits::AU*kUnits::AU)/(radius*radius); // account for 1/r^2
  flux = flux*area*angularView;

  return flux;
  //  return 1*pow( (radius/(1*kUnits::solarRadii)), -2);// Simple r^-2 law
}



// This program takes the detector's position (radius, inclination)
// and returns a double holding the cross section of a cosmic ray 
// background interaction. We want to code softly enough that the 
// function can run in double-pulse and single pulse modes. 


// approximating from https://www.researchgate.net/publication/235341190_SUPERCONDUCTING_TECHNOLOGIES_FOR_THE_ADVANCED_ACCELERATOR_APPLIATIONS_PROGRAM
// we have that protons below 30 MeV will be stopped by 0.1 cm of tunsten


// this program returns the GCR flux at an energy from
// https://alteaspace.wordpress.com/2011/11/27/galactic-cosmic-rays-gcr/
// which gives approximate GCR flux.

double gcrProtonFlux (double radius, double inclination, double energy){
  double intermediate; //holds the value for return
  if(energy < 30*kUnits::MeV){
    intermediate = 0;
  } else if(energy < 40*kUnits::MeV){
    intermediate  = 0.2;
  } else if(energy < 400*kUnits::MeV){
    intermediate  = 0.155*exp(0.006396*energy/kUnits::MeV);
  } else if(energy < 4000*kUnits::MeV){
    intermediate = 2.754*exp(-0.0007998*energy/kUnits::MeV);
  } else if(energy < 10*kUnits::GeV){
    intermediate = 0.3184*exp(-8.066*energy/kUnits::MeV);
  } else {
    intermediate = 0.0001137*exp(-1.2845*energy/kUnits::MeV);
  }
  
  double solidAngle = 4*M_PI;// rougly a solid sphere
  intermediate = intermediate*solidAngle;// answer is in units (m^2 s MeV)^-1
  return intermediate;
}
// total flux
double gcrTotalProtonFlux (double radius, double inclination){
  double Emin = 10*kUnits::MeV;// set to 1/3 energy for protons to pierce shield
  double Emax = 1e7*kUnits::MeV;
  double deltaE = 0.75*kUnits::MeV;
  double integral = 0;
  for (double E = Emin; E<Emax; (E=E+deltaE) ){
    // cout << "I'm on energy value "<< E << "of " << Emax<< endl;
    integral = integral + gcrProtonFlux (radius, inclination, E)*deltaE;
    deltaE = deltaE*sqrt(1.5);
  }



  return integral;
}


// This function produces an approximate cross section of the detector as 
// a function of energy. sqrt(s) =(p1+p2)~ proton momentum
// Assume the oil is ~ 0% hydrogen and 100% carbon by number
// carbon comes from https://www.researchgate.net/publication/48168052_Conceptual_challenges_and_computational_progress_in_X-ray_simulation/figures?lo=1
// hydrogen comes from


double protonCrossSection (double radius, double inclination, double energy){
  double intermediate = 0;
  /*if (energy<0.9*kUnits::MeV){ // Carbon peak
    intermediate = 0.3*0.8496*exp(16.3*energy/kUnits::MeV)*kUnits::barn;
  }
  else if (energy <10*kUnits::MeV){// Carbon change
    intermediate = 0.3*2511488*exp(-0.2530*energy/kUnits::MeV)*kUnits::barn;
  }
  else if (energy <100*kUnits::MeV){// Carbon change
    intermediate = 0.3*239163*exp(-0.01788*energy/kUnits::MeV)*kUnits::barn;
  }
  else{// carbon tail
    0.3*43496*exp(-0.0008378*energy/kUnits::MeV)*kUnits::barn;
  }
  intermediate = intermediate*7/3; //no hydrogen for now
  */
  intermediate = 1*kUnits::barn;// approximation pending investigation
  return intermediate;
}

// This returns the approximate total cross section of the same
double protonTotalCrossSection (double radius, double inclination){
  double integral = 0;
  double Emin = 0.1*kUnits::MeV;
  double Emax = 1e2*kUnits::MeV;
  double deltaE = 0.05*kUnits::MeV; // start with small value
  for (double E = Emin; E<Emax; (E=E+deltaE) ){
    // cout << "I'm on energy value "<< E << "of " << Emax<< endl;
    integral = integral + protonCrossSection (radius, inclination, E)*deltaE;
    deltaE = deltaE*sqrt(1.5);
  }
  
  
  return integral;
}




double cosmicBackground (double radius, double inclination){
  double area = 0.42*M_PI;// meters squared
  /*double Emin = 0.4*kUnits::keV;// set to minimum energy for gallium neutrino interactions
  double Emax = 100*kUnits::MeV;
  Int_t nSteps = 10000;// presently no energy dependence. silly to integrate slowly.
  double deltaE = (Emax-Emin)/nSteps;
  double integral = 0;
  double protonIntegral = 0;
  
  for(double E = Emin; E<Emax; (E=E+deltaE)){
    // Integrates (particleFlux * particleCrossSection  * dE)
    protonIntegral = protonIntegral + gcrProtonFlux(radius, inclination, E)*protonTotalCrossSection(radius, inclination)*deltaE
      + gcrTotalProtonFlux(radius, inclination)*protonCrossSection(radius, inclination, E)*deltaE;   
  }
  // multiplies by the number of targets we have.
  protonIntegral = protonIntegral * protonTargetNumber;// scale to target number
  */
  //integral =  protonIntegral;
  
  double background = area * 3.3*1e4*1.2; // (events/s m^2) * m^2// 120% of the approximate proton rate to account for alphas 
    //integral; // Goes to 0 at 0, and 1 far from the sun
  
  //cout << "cosmicBackground has been called and is returning " << background << " events per second.\n";
  return background; 
}




// This program returns a probability of a radiological decay within a given timespan
double radiologicalBackground (){
  return 123456789; // Goes to 0 at 0, and 1 far from the sun
}

// BACKGROUNDS END -----------------------------------------------------


//
















// This is the main program


// too easy to make a data race in this code
//int nThreads = 12;
//double deltaTimeLimit = timeLimit/nThreads;


int main(){
  User::quadFit qf;
  // These are values that the main program manipulates
  gRandom -> SetSeed(0);
  /*// Moved these so that I can do for each without a data race over the thresholds
    double time, radius, inclination,nuAccumulator,cosmicAccumulator,solarAccumulator=1;
    double nuThresh = (gRandom->Uniform());
    double cosmicThresh = (gRandom->Uniform());
    double solarThresh = (gRandom->Uniform());
  */
  bool isNeutrino, isCosmic, isSolar, isRadio = false; 
  const double cosBack = cosmicBackground(1,1); // currently flat, so we save time
  cout << "cosmic background is " << cosBack << endl;
  string nameOfOutfile = "outfile.root";


  // Control Variables!
  double timeStep = 1*kUnits::s;// time step for iterating
  string myFile = "solar_distance_C3_60(days).csv"; // String that holds r vs time
  double timeLimit = 1*kUnits::day;//1829*kUnits::day; // constant cap
  double detectorRadius = 0.3*kUnits::meter; // approximate value
  double cosmicAcceptance = 0.00001; // 5 nines
  double solarAcceptance = cosmicAcceptance/10000;// four nines better than cosmic
  //double timeLimit = number_of_lines(myFile)*kUnits::day; // soft cap

  // These are the histograms we output

  TH1D* neutrinoRadiusHistogram = new TH1D("neutrinoRadiusHistogram","Neutrino Event Radius",150,3,220);// should always be larger than 0, this program should never have to deal with R>~1.1AU~=230 Rsun
  TH1D* sub35NeutrinoRadiusHistogram = new TH1D("sub35NeutrinoRadiusHistogram","Neutrino Event Radius",150,3,38);// should always be larger than 0, this program should never have to deal with R>~1.1AU~=230 Rsun

  TH1D* neutrinoOnlyRadiusHistogram = new TH1D("neutrinoOnlyRadiusHistogram","Neutrino Only Event Radius",150,3,220);// should always be larger than 0, this program should never have to deal with R>~1.1AU~=230 Rsun
  TH1D* sub35NeutrinoOnlyRadiusHistogram = new TH1D("sub35NeutrinoOnlyRadiusHistogram","Neutrino Only Event Radius",150,3,38);// should always be larger than 0, this program should never have to deal with R>~1.1AU~=230 Rsun
  
  TH1D* cosmicBackRadiusHistogram = new TH1D("cosmicBackRadiusHistogram","Cosmic Background Event Radius",150,3,220);// should always be larger than 0, this program should never have to deal with R>~1.1AU~=230 Rsun
  TH1D* sub35CosmicBackRadiusHistogram = new TH1D("sub35CosmicBackRadiusHistogram","Cosmic Background Event Radius",150,3,38);// should always be larger than 0, this program should never have to deal with R>~1.1AU~=230 Rsun

  TH1D* solarBackRadiusHistogram = new TH1D("solarBackRadiusHistogram","Solar Background Event Radius",150,3,220);// should always be larger than 0, this program should never have to deal with R>~1.1AU~=230 Rsun
  TH1D* sub35SolarBackRadiusHistogram = new TH1D("sub35SolarBackRadiusHistogram","Solar Background Event Radius",150,3,38);// should always be larger than 0, this program should never have to deal with R>~1.1AU~=230 Rsun

  TH1D* radiusHistogram = new TH1D("radiusHistogram","Radius",150,3,220);// always filled no matter what happens

  
  /*
  double* quadraticConstants = johnTheInterpolator(myFile);
  cout << "johnTheInterpolator has spoken!\n\n";
    for (int i = 0;i<5;i++){
      cout << "Day " << (i+1) << " to " << (i+2) << " has constants:\n";
      cout << "c = " << *(quadraticConstants+3*i);
      cout << "b = " << *(quadraticConstants+3*i+1);
      cout << "a = " << *(quadraticConstants+3*i+2) << "\n\n";
    }
    */
 
  ///*
  double closest = 5*kUnits::solarRadii;
  double furthest = 1.1*kUnits::mercurySunDistance;
  double a = (closest+furthest)/2;// semi-major axis
  int nLoops = 100;
  nameOfOutfile = "";
  nameOfOutfile += std::to_string(nLoops);
  nameOfOutfile += "_Orbits.root";
  
  //timeLimit = floor ( ( 2*M_PI*sqrt(a*a*a/( kUnits::G * kUnits::mSun ) ) ) / kUnits::day ) * kUnits::day ;//2*M_PI*sqrt( pow(a,3) / ( kUnits::G*kUnits::mSun ) );
  timeLimit = ( 2*M_PI*sqrt(a*a*a/( kUnits::G * kUnits::mSun ) ) );//2*M_PI*sqrt( pow(a,3) / ( kUnits::G*kUnits::mSun ) );
  
  double* quadraticConstants = qf.elliptical(closest,furthest);
  cout << "elliptical  is done.\n\n";
  for (int i = 0;i<100;i++){
    double c = quadraticConstants[3*i];//*(quadraticConstants+3*i);
    double b = quadraticConstants[3*i];// *(quadraticConstants+3*i+1);
    double a = quadraticConstants[3*i];// *(quadraticConstants+3*i+2);
    //cout << "Day " << (i+1) << " to " << (i+2) << " has constants:\n";
    //cout << "c = " << c;
    //cout << "; b = " << b;
    //cout << "; a = " << a << ".";
    //cout << "Ergo the radius on day " << i+1 << " is " << ( (a*(i+1)*(i+1))+(b*(i+1))+c) << ".\n\n";
  }
  //*/


  //TNtuple* nuNtuple = new TNtuple("nuNtuple","something","Solar_Radii:neutrinoCount:cosmicCount:solarCount");
    double nuAccumulator,cosmicAccumulator,solarAccumulator=0;
    double nuThresh = (gRandom->Uniform());
    double cosmicThresh = (gRandom->Uniform());
    double solarThresh = (gRandom->Uniform());
  for(int i = 0; i < nLoops; i++ ){
    double time, radius, inclination=0;
    time = 0;
    long long int sanityInt = 0;
    cout << "\n\nThis is iteration " << i << ".\n";
    do{
      sanityInt += 1;
      bool doCout = (0 == (sanityInt%20000) );
      

      /* // CSV is slow
      ofstream sanityCSV;
      string sanityName = "sanityCheck";
      sanityName += std::to_string(i);
      sanityName += ".txt";
      sanityCSV.open (sanityName, std::ios::app);
      cout << "The file name is \"" << sanityName << "\".\n";
      */
      if(doCout){
	std::cout << "Before flightPath, time is " << (time/kUnits::days) << " d, radius is " << (radius/kUnits::solarRadii) << " RSol, and inclination is " << inclination <<endl;
      }


      // flightPath(time,radius,inclination);
      radius = (qf.positionMinutes( time, quadraticConstants, timeLimit ));
      
      /*
      cout << "The file name is \"" << sanityName << "\".\n";
      std::cout << "Before flightPath, time is " << (time/kUnits::days) << " d, radius is " << (radius/kUnits::solarRadii) << " RSol, and inclination is " << inclination <<endl;
      sanityCSV << "This is step " << sanityInt << ".\n";
      sanityCSV << radius << ',' << (time/kUnits::hr) << endl ;
      */ // writing to csv is slow

      if (radius > 250*kUnits::solarRadii){
	std::cout << "Before flightPath, time is " << (time/kUnits::days) << " d, radius is " << (radius/kUnits::solarRadii) << " RSol, and inclination is " << inclination <<endl;
	cout << "Oh no! The radius is " << radius/kUnits::solarRadii << " solar radii, which is much further than Earth. I'm setting it to 500 R_Sol or breaking.\n";
	radius = 500*kUnits::solarRadii;
	break;
      }
      if (radius < 1*kUnits::solarRadii) {
	std::cout << "Before flightPath, time is " << (time/kUnits::days) << " d, radius is " << (radius/kUnits::solarRadii) << " RSol, and inclination is " << inclination <<endl;
	cout << "Well shoot, we're at " << radius << " solar radii, which is inside the sun. (Or negative. That makes even less sense). I'm setting it to 0.5 R_sol or breaking.\n";
	  radius = 0.5*kUnits::solarRadii;
	  break;
      }
      
      time = time + timeStep;
      // std::cout << "After flightPath, time is " << time << " radius is " << radius << ", and inclination is " << inclination <<endl;
      
      // The GCR is too big for most timescales directly, but small enough that we don't need it to be
      // I can take the cosmic accumulator, and make it the nearest integer. Once I have that, I can 
      // randomly assign events to bins in the timescale range. A simpler approximation would be to just 
      // assign them to bin 0 through n until I'm out of events. If I do them randomly, then I'll have 
      // to look at the probability of having two hit at once, and the veto performance of that. For
      // the time being, I'll probably stick to even bins, and randomly select a bin for the neutrino.
      bool cosmicBack = true;
      if(cosmicBack){   
	double background = cosmicBackground(radius,1); // checking background
	double vetoedBack = background * cosmicAcceptance * cosmicAcceptance;
	double a = timeStep * cosmicBackground(radius,1) * cosmicAcceptance * cosmicAcceptance;//cosmicBackground(radius,1);
	cosmicAccumulator = cosmicAccumulator + a;

	if(doCout){
	  cout << "Total cosmic backround is " << background << " events per second.\n";
	  cout << "Vetoed cosmic backround is " << vetoedBack << " events per second.\n";
	  std::cout << "The Cosmic Background of " << (radius/kUnits::solarRadii)  << " solarRadii from the center of the sun is : " << a << endl;
	  std::cout << "The Cosmic Background Signal accumulation is " << cosmicAccumulator  << ". " << endl;
	}

	isCosmic = cosmicAccumulator>cosmicThresh;
	if(isCosmic){
	  do{
	    cosmicAccumulator = cosmicAccumulator-cosmicThresh;//0;//cosmicAccumulator - 1;
	    cosmicThresh = (gRandom->Uniform());
	    cosmicBackRadiusHistogram -> Fill( (radius / kUnits::solarRadii) );
	    if(radius < 35*kUnits::solarRadii) sub35CosmicBackRadiusHistogram -> Fill( (radius / kUnits::solarRadii) );
	    //nuNtuple->Fill((radius/kUnits::solarRadii),0,1,0);
	  }
	  while(cosmicAccumulator>cosmicThresh);
	}
      }
      
      bool solarBack = false;
      if(solarBack){  
	double background = solarBackground(radius,1,2*detectorRadius); // checking background
	double vetoedBack = background * solarAcceptance * solarAcceptance;
	double b = timeStep * solarBackground(radius,1,2*detectorRadius)* solarAcceptance * solarAcceptance;
	solarAccumulator = solarAccumulator + b;

	if(doCout){
	  cout << "Total solar backround is " << background << " events per second.\n";
	  cout << "Vetoed solar backround is " << vetoedBack << " events per second.\n";
	  std::cout << "The Solar Background of " << (radius/kUnits::solarRadii)  << " solarRadii from the center of the sun is : " << b << endl;
	  std::cout << "The Solar Background Signal accumulation is " << solarAccumulator  << ". " << endl;
	}

	isSolar = solarAccumulator>solarThresh;
	if(isSolar){
	  if(solarAccumulator>100000) solarAccumulator = 100000;
	  do{
	    long long int weight = 1;
	    double theLog = log10(solarAccumulator);
	    if( theLog > 2){
	      weight = pow(10, floor(theLog - 1));
	      //cout << "I just changed the solar weighting to " << weight << " becase the accumulator was at " << solarAccumulator << endl;
	    }
	    solarAccumulator = solarAccumulator-weight*solarThresh;//0;//solarAccumulator - 1;
	    solarThresh = (gRandom->Uniform());
	    solarBackRadiusHistogram -> Fill( (radius / kUnits::solarRadii), weight );
	    if(radius < 35*kUnits::solarRadii) sub35SolarBackRadiusHistogram -> Fill( (radius / kUnits::solarRadii), weight );
	    //nuNtuple->Fill((radius/kUnits::solarRadii),0,0,1);
	  }
	  while(solarAccumulator>solarThresh);
	}
      }
      
      bool neutrinos = true;
      if(neutrinos){  
	double c = timeStep*neutrinoSignal(radius,1,0.1,20);
	nuAccumulator = nuAccumulator + c;// this creates a data race. I should make a raceless version that is parallel. It may be fast enough to make parallel faster. Each RNG call is ~5ns

	if(doCout){
	  std::cout << "The Neutrino Signal of " << (radius/kUnits::solarRadii)  << " solarRadii from the center of the sun is : " << c << endl;
	  std::cout << "The Neutrino Signal accumulation is " << nuAccumulator  << ". " << endl;
	  std::cout<<endl;
	}

	isNeutrino = nuAccumulator>nuThresh;
	if(isNeutrino){
	  do{
	    bool neutrinoOnly = isNeutrino && !( isCosmic || isSolar || isRadio);
	    nuAccumulator = nuAccumulator-nuThresh;//nuAccumulator-1;
	    nuThresh = (gRandom->Uniform());
	    neutrinoRadiusHistogram -> Fill( (radius / kUnits::solarRadii) );
	    if(radius < 35*kUnits::solarRadii) sub35NeutrinoRadiusHistogram -> Fill( (radius / kUnits::solarRadii) );
	    if(neutrinoOnly) {
	      neutrinoOnlyRadiusHistogram -> Fill( (radius / kUnits::solarRadii) );
	      if(radius < 35*kUnits::solarRadii) sub35NeutrinoOnlyRadiusHistogram -> Fill( (radius / kUnits::solarRadii) );
	    }
	    //nuNtuple->Fill((radius/kUnits::solarRadii),1);
	  }
	  while(nuAccumulator>nuThresh);
	}
      }
      radiusHistogram -> Fill( (radius / kUnits::solarRadii) );
      //if (radius < 15.5*kUnits::solarRadii) break;
    }
    while(time<timeLimit);
  }

  cout << "\n\nThe total orbital period  was " << timeLimit/kUnits::days << " days.\n\n";
  
  // print out the NTuple
  cout << "I'm writing the histograms now.\n";
  const char* nameHolder = nameOfOutfile.c_str();
  TFile *outfile = TFile::Open(nameHolder,"recreate");
  neutrinoRadiusHistogram -> Write();
  sub35NeutrinoRadiusHistogram -> Write();
  neutrinoOnlyRadiusHistogram -> Write();
  sub35NeutrinoOnlyRadiusHistogram -> Write();
  solarBackRadiusHistogram -> Write();
  sub35SolarBackRadiusHistogram -> Write();
  cosmicBackRadiusHistogram -> Write();
  sub35CosmicBackRadiusHistogram -> Write();
  radiusHistogram -> Write();
  outfile -> Close();

  return 0;
}


