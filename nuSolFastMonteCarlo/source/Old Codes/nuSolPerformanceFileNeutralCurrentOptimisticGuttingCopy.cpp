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
#include<math.h>
#include<string>
#include<iostream>
#include<fstream>
#include <iterator>
#include <vector>

#include<TROOT.h>
#include<TFile.h>
#include<TH1D.h>
#include<TRandom.h>
//#include"neutrino.h"
/*
double randomThreshold()
{
   double lower_bound = 0;
   double upper_bound = 1;
   std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
   std::default_random_engine re;
   double a_random_double = unif(re);

   return a_random_double;
   }
*/
using namespace std;




namespace kUnits
{
  // Length (area, volume) Units.
  const double meter = 1;
  const double m = meter;

  const double centimeter = meter/100;
  const double cm = meter/100;

  const double millimeter = meter/1000;
  const double mm = meter/1000;

  const double kilometer = 1000*meter;
  const double km = 1000*meter;

  const double solarRadii = 695700*km;

  const double venusSunDistance = 1.075e8*km;

  const double mercurySunDistance = 4.602e7*km;

  const double AU = 149597870700*m;





  const double barn = 1e-28*m*m;
  const double barns = barn;
  




  const double cc = cm*cm*cm;
  const double milliliter = cc;
  const double mL = cc;

  const double liter = 1000*cc;
  const double L = liter;






  // Time Units
  const double second = 1;
  const double sec = 1;
  const double seconds = 1;
  const double s = 1;

  const double millisecond = s/1000;
  const double milliseconds = s/1000;
  const double ms = s/1000;

  const double microsecond = ms/1000;
  const double microseconds = ms/1000;
  const double us = ms/1000;

  const double nanosecond = us/1000;
  const double nanoseconds = us/1000;
  const double ns = nanosecond;

  const double minute = 60*s;
  const double minutes = 60*s;
  const double min = 60*s;

  const double hour = 60*minute;
  const double hours = 60*minute;
  const double hr = 60*minute;

  const double day = 24*hr;
  const double days = 24*hr;
  const double d = 24*hr;

  const double week = 7*day;
  const double weeks = week;
  const double wk = week;

  const double month = 30*day;
  const double months = month;

  const double year = 365.24219*day;
  const double years = 365.24219*day;
  const double yr = 365.24219*day;
  

  // Energy Units
  const double MeV = 1;

  const double keV = MeV/1;
  
  const double eV = keV/1;
  const double electronVolt = eV;

  const double GeV = 1000*MeV;

  const double Joule = 6.241509e18*eV;
  const double J = Joule;



  // Mass Units
  const double kg = 1;
  const double kilogram = kg;
  
  const double mSun = 1.989e30*kg; // might want to swap with kg for 1
  
  const double g = kg/1000;
  const double gram = g;
  
  const double AMU = kg*1.66053904e-27;

  const double mElectron =  0.51099895000*MeV; // mass of electron (error of 15)

  const double mProton = 938.27208816*MeV; // mass of proton (error of 29)



  // Angle Units
  const double radian = 1;
  const double rad = radian;
  const double radians = radian;

  const double degrees = radian/(2*M_PI )*360;
  const double degree = degrees;
  const double deg = degrees;

  const double arcminutes = degrees/60;
  const double arcminute = arcminutes;
  const double arcmin = arcminutes;

  const double arcseconds = arcminutes/60;
  const double arcsecond = arcseconds;
  const double arcsec = arcseconds;
  
  // Physical Constants
  const double c = 299792458*m/s;
  const double speedOfLight = c;

  const double G = 6.67430e-11; // this one's units are all 1 in this system
  const double newtonGravity = G; 

  const double reducedFermiConstant = 1.1663787e-5/GeV/GeV;
  const double G_F = reducedFermiConstant;

  const double weakCoupling = 0.22290; // sin^2 of weak mixing angle; error of 30 in constant
  const double theta_W = asin(sqrt(0.22290)); // Weak Mixing/Weinberg angle; error of 30 in constant




}

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
    //  cout << "Line " << count << " gives datum " << datum <<endl;
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
      /*( 2*y[i-1] - y[i] - y[i+1] ) / ( -3*i ) 
	- a * ( -6*i-1 ) / ( -3 );*/
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

double positionDays(double time, double* quadConst, double maxTime){
  //cout << "Position function has been called.\n";
  double currentRadius=0;
  int day = floor(time/kUnits::day);
  /* cout is slow
    cout << "The time is " << (time/kUnits::hr) << "hours in, and we are " 
       << day << " days in.\n";
  */
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
    /* couts are slow
    cout << "a = " << *( quadConst+2  + 3*( day-1 ) ) << ". b = "
	 << *( quadConst+1  + 3*( day-1 ) ) << ". c = "
	 << *(quadConst + 3*( day-1 ) ) << ".\n";
    */

    currentRadius = *(quadConst + 3*( day-1 ) )
      + *( quadConst+1  + 3*( day-1 ) ) * time /kUnits::hr 
      + *( quadConst+2  + 3*( day-1 ) ) * time * time /kUnits::hr /kUnits::hr;
    
  }
  //cout << "The current radius is " << currentRadius << endl;
  return currentRadius*kUnits::solarRadii; // return the radius in solar radii
}

// minutes version
double positionMinutes(double time, double* quadConst, double maxTime){
  cout << "Position function has been called.\n";
  double currentRadius=0;
  int day = floor(time/kUnits::day);
  //cout << "The time is " << (time/kUnits::hr) << "hours in, and we are " 
  //     << day << " days in.\n";
  if(time < (1*kUnits::min)){
    currentRadius = *(quadConst)+time*(*( quadConst+1 ) )/kUnits::min 
      + time * time * (*( quadConst+2 ) ) /kUnits::min /kUnits::min;
  }
  else if (time > (maxTime - kUnits::min)){ // last day from earlier data
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
    /* couts are slow
      cout << "a = " << *( quadConst+2  + 3*( day-1 ) ) << ". b = "
           << *( quadConst+1  + 3*( day-1 ) ) << ". c = "
	   << *(quadConst + 3*( day-1 ) ) << ".\n";
    */
    currentRadius = *(quadConst + 3*( day-1 ) )
      + *( quadConst+1  + 3*( day-1 ) ) * time /kUnits::min 
      + *( quadConst+2  + 3*( day-1 ) ) * time * time /kUnits::min /kUnits::min;
    
  }
  //cout << "The current radius is " << currentRadius << endl;
  return currentRadius*kUnits::solarRadii; // return the radius in solar radii
}


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
  double energy = slope*radius + intercept;
  //cout << "The energy is " << energy/kUnits::MeV << " Mev at radius " << radius/kUnits::AU << " AU, the intercept is " << intercept/kUnits::MeV << "MeV and the slope is " << slope/kUnits::MeV*kUnits::AU << " MeV/AU.\n";
  return energy;
}

// NEUTRINOS BY PROCESS AND ENERGY


// from https://arxiv.org/pdf/0804.3899.pdf
double diffCrossSec(){
  double naturalUnits = (1 // high energy left-left interaction
			 + (1 - 2*kUnits::weakCoupling + 4/3*kUnits::weakCoupling*kUnits::weakCoupling) // Charge Current right-left
			 +(1/4 - kUnits::weakCoupling + 4/3*kUnits::weakCoupling*kUnits::weakCoupling) // Neutral Current right-left
			 )*kUnits::G_F*kUnits::G_F*2*kUnits::mElectron/M_PI;
  double conversion = 1/(5.06*5.06e26)*kUnits::cm*kUnits::cm;// conversion from GeV^-2 to cm^2
  double mksUnits = naturalUnits*conversion;
  return mksUnits;
}



double hepRate (double eMin, double eMax){
  eMin = eMin/kUnits::MeV;// fix units for integrating
  eMax = eMax/kUnits::MeV;
  
  // 5 regions with 5 fits
  double gMin = 0.01; // hard to detect <10 keV
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
  integral += const1 * ( pow(pow1Max,pow1+2) - pow(pow1Min,pow1+2) ) 
              / (pow1+2) * (1 + 1 / (pow1+1) ); // Power integral
  integral += const2 * ( pow(pow2Max,pow2+2) - pow(pow2Min,pow2+2) ) / (pow2+2) * (1 + 1 / (pow2+1) );
  integral += const3 * ( pow(pow3Max,pow3+2) - pow(pow3Min,pow3+2) ) / (pow3+2) * (1 + 1 / (pow3+1) );
  integral += constant * (constMax*constMax-constMin*constMin);// constant integral
  
  integral = integral*diffCrossSec();

  return integral/(kUnits::cm)/(kUnits::cm)/(kUnits::s)/(M_PI*4);// per second per cm^2 per steradian
}

double B8Rate (double eMin, double eMax){
  eMin = eMin/kUnits::MeV;// fix units for integrating
  eMax = eMax/kUnits::MeV;
  
  // 5 regions with 5 fits
  double gMin = 0.01; // hard to detect <10 keV
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
  integral += const1 * ( pow(pow1Max,pow1+2) - pow(pow1Min,pow1+2) ) 
              / (pow1+2) * (1 + 1 / (pow1+1) ); // Power integral
  integral += const2 * ( pow(pow2Max,pow2+2) - pow(pow2Min,pow2+2) ) / (pow2+2) * (1 + 1 / (pow2+1) );
  integral += const3 * ( pow(pow3Max,pow3+2) - pow(pow3Min,pow3+2) ) / (pow3+2) * (1 + 1 / (pow3+1) );
  integral += constant * (constMax*constMax-constMin*constMin);// constant integral
  
  integral = integral*diffCrossSec();
  
  return integral/(kUnits::cm)/(kUnits::cm)/(kUnits::s)/(M_PI*4);// per second per cm^2 per steradian
}

double F17Rate (double eMin, double eMax){
  eMin = eMin/kUnits::MeV;// fix units for integrating
  eMax = eMax/kUnits::MeV;
  
  // 5 regions with 5 fits
  double gMin = 0.01; // hard to detect <10 keV
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
  integral += const1 * ( pow(pow1Max,pow1+2) - pow(pow1Min,pow1+2) ) 
              / (pow1+2) * (1 + 1 / (pow1+1) ); // Power integral
  integral += const2 * ( pow(pow2Max,pow2+2) - pow(pow2Min,pow2+2) ) / (pow2+2) * (1 + 1 / (pow2+1) );
  integral += const3 * ( pow(pow3Max,pow3+2) - pow(pow3Min,pow3+2) ) / (pow3+2) * (1 + 1 / (pow3+1) );
  integral += constant * (constMax*constMax-constMin*constMin);// constant integral
  
  integral = integral*diffCrossSec();
  
  return integral/(kUnits::cm)/(kUnits::cm)/(kUnits::s)/(M_PI*4);// per second per cm^2 per steradian
}

double O15Rate (double eMin, double eMax){
  eMin = eMin/kUnits::MeV;// fix units for integrating
  eMax = eMax/kUnits::MeV;
  
  // 5 regions with 5 fits
  double gMin = 0.01; // hard to detect <10 keV
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
  integral += const1 * ( pow(pow1Max,pow1+2) - pow(pow1Min,pow1+2) ) 
              / (pow1+2) * (1 + 1 / (pow1+1) ); // Power integral
  integral += const2 * ( pow(pow2Max,pow2+2) - pow(pow2Min,pow2+2) ) / (pow2+2) * (1 + 1 / (pow2+1) );
  integral += const3 * ( pow(pow3Max,pow3+2) - pow(pow3Min,pow3+2) ) / (pow3+2) * (1 + 1 / (pow3+1) );
  integral += constant * (constMax*constMax-constMin*constMin);// constant integral
  
  integral = integral*diffCrossSec();
  
  return integral/(kUnits::cm)/(kUnits::cm)/(kUnits::s)/(M_PI*4);// per second per cm^2 per steradian
}

double N13Rate (double eMin, double eMax){
  eMin = eMin/kUnits::MeV;// fix units for integrating
  eMax = eMax/kUnits::MeV;
  
  // 5 regions with 5 fits
  double gMin = 0.01; // hard to detect <10 keV
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
  integral += const1 * ( pow(pow1Max,pow1+2) - pow(pow1Min,pow1+2) ) 
              / (pow1+2) * (1 + 1 / (pow1+1) ); // Power integral
  integral += const2 * ( pow(pow2Max,pow2+2) - pow(pow2Min,pow2+2) ) / (pow2+2) * (1 + 1 / (pow2+1) );
  integral += const3 * ( pow(pow3Max,pow3+2) - pow(pow3Min,pow3+2) ) / (pow3+2) * (1 + 1 / (pow3+1) );
  integral += constant * (constMax*constMax-constMin*constMin);// constant integral
  
  integral = integral*diffCrossSec();
  
  return integral/(kUnits::cm)/(kUnits::cm)/(kUnits::s)/(M_PI*4);// per second per cm^2 per steradian
}

double ppRate (double eMin, double eMax){
  eMin = eMin/kUnits::MeV;// fix units for integrating
  eMax = eMax/kUnits::MeV;
  
  // 5 regions with 5 fits
  double gMin = 0.01; // hard to detect <10 keV
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
  integral += const1 * ( pow(pow1Max,pow1+2) - pow(pow1Min,pow1+2) ) 
              / (pow1+2) * (1 + 1 / (pow1+1) ); // Power integral
  integral += const2 * ( pow(pow2Max,pow2+2) - pow(pow2Min,pow2+2) ) / (pow2+2) * (1 + 1 / (pow2+1) );
  integral += const3 * ( pow(pow3Max,pow3+2) - pow(pow3Min,pow3+2) ) / (pow3+2) * (1 + 1 / (pow3+1) );
  integral += constant * (constMax*constMax-constMin*constMin);// constant integral
  
  integral = integral*diffCrossSec();
  
  return integral/(kUnits::cm)/(kUnits::cm)/(kUnits::s)/(M_PI*4);// per second per cm^2 per steradian
}

double pepRate(double eMin, double eMax){
  double eDelta = 1.3*kUnits::MeV; //placeholder
  if ( (eMin < eDelta) && (eMax > eDelta)){
    return 1.3e8/(kUnits::cm)/(kUnits::cm)/(kUnits::s)/(M_PI*4)*diffCrossSec()*eDelta;// per second per cm^2 per steradian
  }
  else{
    return 0;
  } 
}

double Be7Rate(double eMin, double eMax){
  double eDelta1 = 0.38*kUnits::MeV; //placeholder
  double eDelta2 = 0.87*kUnits::MeV; //placeholder
  if ( (eMin < eDelta1) && (eMax > eDelta1)){
    return 5e8/(kUnits::cm)/(kUnits::cm)/(kUnits::s)/(M_PI*4)*diffCrossSec()*eDelta1;// per second per cm^2 per steradian
  }
  else if ( (eMin < eDelta2) && (eMax > eDelta2)){
    return 4.1e9/(kUnits::cm)/(kUnits::cm)/(kUnits::s)/(M_PI*4)*diffCrossSec()*eDelta2;// per second per cm^2 per steradian
  }
  else{
    return 0;
  } 
}


double detectorDensity = 1*kUnits::g/pow(kUnits::cm,3); // plastic density
double kgDetector = 1000;
double detectorVolume = kgDetector/detectorDensity;

double nProtons = 6+2+6+1+6*6+4+6+3 ; // based on Polyvinyl Toluene (base for a scintillator) (CH2CH(C6H4CH3))
double massMolecules = (12+2+12+1+6*12+4+12+3)*kUnits::AMU;
double electronNumber = kgDetector*nProtons/massMolecules; // det mass * protons/mass

// This program takes the detector's position (radius, inclination)
// and returns a double holding the flux of  neutrinos at that
// point. We want to code softly enough that the function
// can run in double-pulse and single pulse modes. 

double neutrinoSignal (double radius, double inclination, double eMin, double eMax){
  
  // number of interactions per second per atom at 1 AU. 
  // Calculated from effective cross sections and solar
  // neutrino flux by source graph. 
  double neutrinoInteractionRate = 0;//1.05957e-34;
  neutrinoInteractionRate += ppRate(eMin, eMax);
  neutrinoInteractionRate += pepRate(eMin, eMax);
  neutrinoInteractionRate += Be7Rate(eMin, eMax);
  neutrinoInteractionRate += N13Rate(eMin, eMax);
  neutrinoInteractionRate += O15Rate(eMin, eMax);
  neutrinoInteractionRate += F17Rate(eMin, eMax);
  neutrinoInteractionRate += B8Rate(eMin, eMax);
  neutrinoInteractionRate += hepRate(eMin, eMax);
 
  

  // find the rate per electron
  neutrinoInteractionRate = neutrinoInteractionRate*electronNumber;
  
  // find the rate at arbitrary radius
  neutrinoInteractionRate = neutrinoInteractionRate*pow( (radius/kUnits::AU),-2);
  
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
  return (1e-5 / pow(kUnits::cm,2))/kUnits::second/pow(radius/kUnits::solarRadii,2);// near zero outsie of solar flares?
}

// This program takes the detector's position (radius, inclination)
// and returns a double holding the cross section of a solar wind 
// background interaction. We want to code softly enough that the 
// function can run in double-pulse and single pulse modes.

double neutronTargetNumber = 1e27; // 10 times gallium target number for now
double protonTargetNumber = 1e28; // ~assuming the target is 400 pounds of tungsten

double detectorRadius = pow(detectorVolume*3/4/M_PI, 1/3); // V = 4/3 pi r^3
double detectorArea = M_PI*detectorRadius*detectorRadius; // A = pi r^2


double solarBackground (double radius, double inclination){
  // assuming the tungsten is 5 cm thick, the stopping power is 5*20
  // multiplied by the data at https://physics.nist.gov/cgi-bin/Star/ap_table.pl
  // This means we only need to worry about protons >~ 10 MeV
  // We'll assume the 10 MeV warning threshold from the GOES data,
  // 10 cm^-2 s^-1 sr^-1, as a basis for analysis.
  // If this is too high, we can look more in depth for solar maximum data.
  // https://www.swpc.noaa.gov/products/goes-proton-flux-dynamic-plot
  //double area = M_PI*(0.3*kUnits::meter)*(0.3*kUnits::meter);// r is about 0.2 meters for our detector
  // but there's also some small bitys of the sides that will be hit

  double angularView = M_PI;// quarter of view, probably alwasys smaller than this, but dependent on radius
                              // this might take the 1/r^2 depndence term in later.
  
  double earthValue = 10/kUnits::cm/kUnits::cm/kUnits::sec;// sr^-1
  double flux = earthValue *(kUnits::AU*kUnits::AU)/(radius*radius); // account for 1/r^2
  flux = flux*detectorArea*angularView;

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
  //double area = 0.42*M_PI;// meters squared
  
  double background = detectorArea * 3.3*1e4*1.2; // (events/s m^2) * m^2// 120% of the approximate proton rate to account for alphas 
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
  // These are values that the main program manipulates
  gRandom -> SetSeed(0);
  bool isNeutrino, isCosmic, isSolar, isRadio = false; 
  const double cosBack = cosmicBackground(1,1); // currently flat, so we save time
  cout << "cosmic background is " << cosBack << endl;
  string nameOfOutfile = "outfile.root";
  int nLoops = 100;


  // Control Variables!
  bool inward = false;
  double timeStep = 1*kUnits::min;// time step for iterating
  string myFile = "../outwardTrajectory.txt";//solar_distance_C3_60(days).csv; // String that holds r vs time
  double timeLimit = 1*kUnits::day;//1829*kUnits::day; // constant cap
  double cosmicAcceptance = 0.00001; // 5 nines
  double solarAcceptance = cosmicAcceptance;// same as cosmic
  //double timeLimit = number_of_lines(myFile)*kUnits::day; // soft cap

  // Derived Control Variables
  double histoLowerLimit = 0;
  if(inward) histoLowerLimit = 3;
  else histoLowerLimit = 210; // lower limit in solar radii
  double histoUpperLimit = 0;
  if(inward) histoUpperLimit = 220 ;
  else histoUpperLimit = 6450;
  double histoBins = 0;
  if(inward) histoBins = 150;
  else histoBins = 1000;

  // These are the histograms we output

  TH1D* neutrinoRadiusHistogram = new TH1D("neutrinoRadiusHistogram","Neutrino Event Radius",histoBins,histoLowerLimit,histoUpperLimit);// should always be larger than 0, this program should never have to deal with R>~1.1AU~=230 Rsun
  TH1D* sub35NeutrinoRadiusHistogram = new TH1D("sub35NeutrinoRadiusHistogram","Neutrino Event Radius",150,3,38);// should always be larger than 0, this program should never have to deal with R>~1.1AU~=230 Rsun

  TH1D* neutrinoOnlyRadiusHistogram = new TH1D("neutrinoOnlyRadiusHistogram","Neutrino Only Event Radius",histoBins,histoLowerLimit,histoUpperLimit);// should always be larger than 0, this program should never have to deal with R>~1.1AU~=230 Rsun
  TH1D* sub35NeutrinoOnlyRadiusHistogram = new TH1D("sub35NeutrinoOnlyRadiusHistogram","Neutrino Only Event Radius",150,3,38);// should always be larger than 0, this program should never have to deal with R>~1.1AU~=230 Rsun
  
  TH1D* cosmicBackRadiusHistogram = new TH1D("cosmicBackRadiusHistogram","Cosmic Background Event Radius",histoBins,histoLowerLimit,histoUpperLimit);// should always be larger than 0, this program should never have to deal with R>~1.1AU~=230 Rsun
  TH1D* sub35CosmicBackRadiusHistogram = new TH1D("sub35CosmicBackRadiusHistogram","Cosmic Background Event Radius",150,3,38);// should always be larger than 0, this program should never have to deal with R>~1.1AU~=230 Rsun

  TH1D* solarBackRadiusHistogram = new TH1D("solarBackRadiusHistogram","Solar Background Event Radius",histoBins,histoLowerLimit,histoUpperLimit);// should always be larger than 0, this program should never have to deal with R>~1.1AU~=230 Rsun
  TH1D* sub35SolarBackRadiusHistogram = new TH1D("sub35SolarBackRadiusHistogram","Solar Background Event Radius",150,3,38);// should always be larger than 0, this program should never have to deal with R>~1.1AU~=230 Rsun

  TH1D* radiusHistogram = new TH1D("radiusHistogram","Radius",histoBins,histoLowerLimit,histoUpperLimit);// always filled no matter what happens

  TH1D* neutrinoFocusHistogram = new TH1D("neutrinoFocusHistogram","Radius (rSol)",histoBins,histoLowerLimit,histoUpperLimit);// always filled no matter what happens
  TH1D* neutrinoOnlyFocusHistogram = new TH1D("neutrinoOnlyFocusHistogram","Radius (AU)",histoBins,histoLowerLimit,histoUpperLimit);// always filled no matter what happens

  
  
  double* quadraticConstants = johnTheInterpolator(myFile);
  cout << "johnTheInterpolator has spoken!\n\n";
    for (int i = 0;i<5;i++){
      cout << "Day " << (i+1) << " to " << (i+2) << " has constants:\n";
      cout << "c = " << *(quadraticConstants+3*i);
      cout << "b = " << *(quadraticConstants+3*i+1);
      cout << "a = " << *(quadraticConstants+3*i+2) << "\n\n";
    }
    timeLimit = kUnits::days*number_of_lines(myFile);
  
  
  


  
  for(int i = 0; i < nLoops; i++ ){
    double time, radius, inclination,nuAccumulator,nuAccumulator2,cosmicAccumulator,solarAccumulator=1;
    time = 0;
    int sanityInt = 0;
    double nuThresh = (gRandom->Uniform());
    double nuThresh2 = (gRandom->Uniform());
    double cosmicThresh = (gRandom->Uniform());
    double solarThresh = (gRandom->Uniform());
    cout << "\n\nThis is iteration " << i << ".\n";
    do{
      sanityInt += 1;
      bool doCout = false; //(0 == (sanityInt%200000) );
      
      if(doCout){
	std::cout << "Before flightPath, time is " << (time/kUnits::days) << " d, radius is " << (radius/kUnits::solarRadii) << " RSol, and inclination is " << inclination <<endl;
      }


      radius = (positionDays (time, quadraticConstants, timeLimit )); // for file
      double beforeRadius = (positionDays ( (time-timeStep/2) , quadraticConstants, timeLimit )); // for file
      double afterRadius = (positionDays ( (time+timeStep/2) , quadraticConstants, timeLimit )); // for file
      
    if (radius > 30*kUnits::AU){
	std::cout << "Before flightPath, time is " << (time/kUnits::days) << " d, radius is " << (radius/kUnits::solarRadii) << " RSol, and inclination is " << inclination <<endl;
	cout << "Oh no! The radius is " << radius/kUnits::AU << " AU, which is stupid far. I'm setting it to 500 R_Sol or breaking.\n";
	radius = 500*kUnits::solarRadii;
	break;
      }
      if (radius < 1) {
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
	double cosmicGeoAccpet = 4*M_PI*kUnits::solarRadii*kUnits::solarRadii/(4/M_PI)/(radius*radius);// fraction of the sky that could be solar. We'll assume that we can isolate to 2 rSol (unrealistically good)
	double background = cosmicBackground(radius,1); // checking background
	double vetoedBack = background * cosmicAcceptance * cosmicGeoAccpet;
	double a = timeStep * cosmicBackground(radius,1) * cosmicAcceptance * cosmicGeoAccpet;
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
      
      bool solarBack = true;
      if(solarBack){  
	double background = solarBackground(radius,1); // checking background
	double vetoedBack = background * solarAcceptance;
	double b = timeStep * solarBackground(radius,1)* solarAcceptance;
	solarAccumulator = solarAccumulator + b;

	if(doCout){
	  cout << "Total solar backround is " << background << " events per second.\n";
	  cout << "Vetoed solar backround is " << vetoedBack << " events per second.\n";
	  std::cout << "The Solar Background of " << (radius/kUnits::solarRadii)  << " solarRadii from the center of the sun is : " << b << endl;
	  std::cout << "The Solar Background Signal accumulation is " << solarAccumulator  << ". " << endl;
	}

	isSolar = solarAccumulator>solarThresh;
	if(isSolar){
	  do{
	  solarAccumulator = solarAccumulator-solarThresh;//0;//solarAccumulator - 1;
	  solarThresh = (gRandom->Uniform());
	  solarBackRadiusHistogram -> Fill( (radius / kUnits::solarRadii) );
	  if(radius < 35*kUnits::solarRadii) sub35SolarBackRadiusHistogram -> Fill( (radius / kUnits::solarRadii) );
	  //nuNtuple->Fill((radius/kUnits::solarRadii),0,0,1);
	  }
	  while(solarAccumulator>solarThresh);
	}
      }
      
      bool neutrinos = true;
      if(neutrinos){  
	double c = timeStep*neutrinoSignal(radius,1,0.2,17);
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
      if(neutrinos){  
	double d = 880*timeStep*neutrinoSignal(radius,1,rToE(beforeRadius),rToE(afterRadius))*pow( (radius/kUnits::AU),2);// 880x and restore the 1/r^2 dividede out for the focus stuff
	nuAccumulator2 = nuAccumulator2 + d;// this creates a data race. I should make a raceless version that is parallel. It may be fast enough to make parallel faster. Each RNG call is ~5ns
	
	if(doCout){
	  std::cout << "The Neutrino Signal at the lensing focus is : " << d << endl;
	  std::cout << "The Neutrino Signal accumulation is " << nuAccumulator2  << ". " << endl;
	  std::cout<<endl;
	}

	isNeutrino = nuAccumulator2>nuThresh2;
	if(isNeutrino){
	  bool neutrinoOnly = isNeutrino && !( isCosmic || isSolar || isRadio);
	  do{
	    nuAccumulator2 = nuAccumulator2-nuThresh2;//nuAccumulator-1;
	    nuThresh2 = (gRandom->Uniform());
	    neutrinoFocusHistogram-> Fill( (radius / kUnits::solarRadii) );
	    if(neutrinoOnly) {
	      neutrinoOnlyFocusHistogram-> Fill( (radius / kUnits::solarRadii) );
	    }
	  }
	  while(nuAccumulator2>1);
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
  neutrinoFocusHistogram -> Write();
  neutrinoOnlyFocusHistogram -> Write();
  solarBackRadiusHistogram -> Write();
  cosmicBackRadiusHistogram -> Write();
  sub35CosmicBackRadiusHistogram -> Write();
  radiusHistogram -> Write();
  outfile -> Close();

  return 0;
}


