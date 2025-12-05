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
#include <iterator>
#include <vector>

#include"TROOT.h"
#include"TRandom.h"
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

  const double venusSunDistance = 1.075*pow(10,8)*km;

  const double mercurySunDistance = 4.602*pow(10,7)*km;

  const double AU = 149597870700*m;





  const double barn = pow(10,-28)*m*m;
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



  // Mass Units
  const double kg = 1;
  const double kilogram = kg;
  
  const double mSun = 1.989*pow(10,30)*kg; // might want to swap with kg for 1



  // Energy Units
  const double MeV = 1;

  const double keV = MeV/1;
  
  const double eV = keV/1;
  const double electronVolt = eV;

  const double GeV = 1000*MeV;

  const double Joule = 6.241509*pow(10,18)*eV;
  const double J = Joule;



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

  const double G = 6.67430*pow(10,-11); // units are all 1 in this system
  const double newtonGravity = G; 




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
  std::ifstream data;
  data.open (theFile);
  cout << "I have opened the file \"" << theFile << "\"\n\n";
  std::string line;
  int nLines=0;
  while (std::getline(data, line))
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
  int cells = ceil(timeMax/kUnits::minute);
  //double theValues[cells];
  double theValues[cells];
  int sanityCount = 0;


  cout << "I'm about to find the position/radii/etc for the ellipse uising " << intSteps << " steps. I am going to put them into " << cells << " cells.\n\n";
  for (int i = 0; i < intSteps; i++){
    int cellNum = ceil(i*timeStep/kUnits::minute);    
    //    cout << "day = " << i*timeStep/kUnits::days << endl;
    if (i%100==0){
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
    if(i%100==0){
      cout << "The array value for the " << i << "th radius is r = " << theValues[i] << endl;
    }
  }


  cout << "\n\nI'm about to find the fit parameters for the ellipse, which should hold " << (3*(cells-1))  << "values.\n\n";
  double* theReturns = new double[(3*(cells-1))];// return the values as one long array
  // this type of declaration shuld fix a memory issue. We keep theReturns around until deleted


  for(int i=1; i<(cells-1);i++){
    // values found from calculating (a,b,c)
    // =inverse( (1,x(i-1),x^2(i-1); 1,x(i),x^2(i)); 1,x(i+1),x^2(i+1) )
    // (inner product) (c,b,a)^T
    double a = theValues[i-1]/2-theValues[i]+theValues[i+1]/2;
    // ( 2*theValues[i]-theValues[i-1]-theValues[i+1] ) // numerator
    // / ( -2 ); // denominator
    double b = ( (-2*i-1)*theValues[i-1] + 4*i*theValues[i] + (-2*i+1)*theValues[i+1] )/2;
    /*( 2*theValues[i-1] - theValues[i] - theValues[i+1] ) / ( -3*i ) 
      - a * ( -6*i-1 ) / ( -3 );*/
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
  cout << "Position function has been called.\n";
  double currentRadius=0;
  int day = floor(time/kUnits::day);
  cout << "The time is " << (time/kUnits::hr) << "hours in, and we are " 
       << day << " days in.\n";
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
  cout << "The current radius is " << currentRadius << endl;
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
    cout << "a = " << *( quadConst+2  + 3*( day-1 ) ) << ". b = "
	 << *( quadConst+1  + 3*( day-1 ) ) << ". c = "
	 << *(quadConst + 3*( day-1 ) ) << ".\n";
    currentRadius = *(quadConst + 3*( day-1 ) )
      + *( quadConst+1  + 3*( day-1 ) ) * time /kUnits::min 
      + *( quadConst+2  + 3*( day-1 ) ) * time * time /kUnits::min /kUnits::min;
    
  }
  cout << "The current radius is " << currentRadius << endl;
  return currentRadius*kUnits::solarRadii; // return the radius in solar radii
}


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
double kgGallium = 100;
double galliumAtomNumber =kgGallium * 8.637*pow(10,24);// 1 kg of gallium for now

// This program takes the energy of an incoming neutrino and outputs
// a cross section for the gallium interaction at that energy.

/*double neutrinoCrossSection (double radius, double inclination, double energy){
  return (pow(10,-46) * 8.6*pow((kUnits::cm),2));// approximately constant for all energies.
}

double neutrinoTotalCrossSection (double radius, double inclination){
  double Emin = 234.7*kUnits::keV;// set to minimum energy for gallium neutrino interactions + 1 sigma (1.2) https://arxiv.org/pdf/1710.06326.pdf
  double Emax = 10*kUnits::MeV;// relatively small currently because we don't know energy dependence
  Int_t nSteps = 1;// presently no energy dependence. silly to integrate slowly.
  double deltaE = (Emax-Emin)/nSteps;
  double integral = 0;
  for(double E = Emin; E<Emax; (E=E+deltaE)){
    integral = integral + neutrinoCrossSection(radius, inclination, E)*deltaE;
  }
  
  return integral;// approximately constant for all energies.
}
    

// This program takes position as a function of radius and inclination,
// and also the energy of a neutrino to find the neutrino flux per unit 
// energy for that energy and position.

double neutrinoFlux (double radius, double inclination, double energy){
  return (pow(10,11))*(1/pow( (radius/(1*kUnits::AU)), 2)); 
    // Simple (r^-2 law * E)
}

double neutrinoTotalFlux (double radius, double inclination){
  double Emin = 234.7*kUnits::keV;// set to minimum energy for gallium neutrino interactions + 1 sigma (1.2) https://arxiv.org/pdf/1710.06326.pdf
  double Emax = 10*kUnits::MeV;// relatively small currently because we don't know energy dependence
  Int_t nSteps = 1;// presently no energy dependence. silly to integrate slowly.
  double deltaE = (Emax-Emin)/nSteps;
  double integral = 0;
  for(double E = Emin; E<Emax; (E=E+deltaE)){
    integral = integral + neutrinoFlux(radius, inclination, E);
  }
  
  return integral;// approximately constant for all energies.
    // Simple (r^-2 law * E)
}
*/


// all fluxes come from the famous specturm downloaded from 
// https://falcon.phy.queensu.ca/SNO+/about/solar-neutrinos.html
// pixel counted in GIMP, and fitted to exponential segments in
// OpenOffice

double ppFlux(double radius, double inclination){
  return 1;
}






// This program takes the detector's position (radius, inclination)
// and returns a double holding the flux of  neutrinos at that
// point. We want to code softly enough that the function
// can run in double-pulse and single pulse modes. 

double neutrinoSignal (double radius, double inclination){
  /*double Emin = 234.7*kUnits::keV;// set to minimum energy for gallium neutrino interactions + 1 sigma (1.2) https://arxiv.org/pdf/1710.06326.pdf
  double Emax = 100*kUnits::MeV;
  Int_t nSteps = 20;// presently no energy dependence. silly to integrate slowly.
  double deltaE = (Emax-Emin)/nSteps;
  double integral = 0;
  
  for(double E = Emin; E<Emax; (E=E+deltaE)){
    
    // Integrates (neturinoFlux * neutrinoCrossSection  * dE)
    integral = integral + neutrinoTotalFlux(radius, inclination)*neutrinoCrossSection(radius, inclination, E)*deltaE
      + neutrinoFlux(radius, inclination, E)*neutrinoTotalCrossSection(radius, inclination)*deltaE;   
  }
  // multiplies by the number of targets we have.
  integral = integral * galliumAtomNumber;
  
  double totalNeutrinoInteractionRate= 0;

  double ppCrossSec = 11.72*pow(10,-46)*kUnits::cm*kUnits::cm;
  double pepCrossSec = 204*pow(10,-46)*kUnits::cm*kUnits::cm;
  double Be7CrossSec = 71.7*pow(10,-46)*kUnits::cm*kUnits::cm;
  double N13CrossSec = 60.4*pow(10,-46)*kUnits::cm*kUnits::cm;
  double O15CrossSec = 113.7*pow(10,-46)*kUnits::cm*kUnits::cm;
  double F17CrossSec = 113.9*pow(10,-46)*kUnits::cm*kUnits::cm;
  double Ar37CrossSec = 70.0*pow(10,-46)*kUnits::cm*kUnits::cm;
  double Cr51CrossSec = 58.1*pow(10,-46)*kUnits::cm*kUnits::cm;

  double B8CrossSec = 2.40*pow(10,-42)*kUnits::cm*kUnits::cm;
  double hepCrossSec = 7.14*pow(10,-42)*kUnits::cm*kUnits::cm;
  
  double Emin = 234.7*kUnits::keV;// set to minimum energy for gallium neutrino interactions + 1 sigma (1.2) https://arxiv.org/pdf/1710.06326.pdf
  double Emax = 20*kUnits::MeV;
  double deltaE = 50*kUnits::keV;
  double ppIntegral = 0;

  for (double E = Emin; E<Emax; (E=E+deltaE) ){
    // cout << "I'm on energy value "<< E << "of " << Emax<< endl;
    ppIntegral = ppIntegral + gcrProtonFlux (radius, inclination, E)*deltaE;
    deltaE = deltaE*sqrt(1.5);
  }
  */


  // number of interactions per second per atom at 1 AU. 
  // Calculated from effective cross sections and solar
  // neutrino flux by source graph. 
  double neutrinoInteractionRate = 1.05957*pow(10,-34);


  // find the rate per kg of gallium
  neutrinoInteractionRate = neutrinoInteractionRate*galliumAtomNumber;
  
  // find the rate at arbitrary radius
  neutrinoInteractionRate = neutrinoInteractionRate*pow( (radius/kUnits::AU),-2);
  
  return neutrinoInteractionRate;
}

// NEUTRINOS END ---------------------------------------------------------








// BACKGROUNDS START ----------------------------------------------------


// This function determines the neutron cross section

double neutronCrossSection (double radius, double inclination, double energy){
  return (pow(10,1) *kUnits::barns);// approximately constant or smaller than this for all energies of interest.
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
  return (pow(10,-5) / pow(kUnits::cm,2))/kUnits::second/pow(radius/kUnits::solarRadii,2);// near zero outsie of solar flares?
}

// This program takes the detector's position (radius, inclination)
// and returns a double holding the cross section of a solar wind 
// background interaction. We want to code softly enough that the 
// function can run in double-pulse and single pulse modes.

double neutronTargetNumber = pow(10,27); // 10 times gallium target number for now
double protonTargetNumber = pow(10,28); // ~assuming the target is 400 pounds of tungsten

double solarBackground (double radius, double inclination){
  double Emin = 0.4*kUnits::keV;// set to minimum energy for gallium neutrino interactions
  double Emax = 100*kUnits::MeV;
  Int_t nSteps = 100;// presently no energy dependence. silly to integrate slowly.
  double deltaE = (Emax-Emin)/nSteps;
  double neutronIntegral = 0;
  
  for(double E = Emin; E<Emin+1/*Emax*/; (E=E+deltaE)){
    // Integrates (particleFlux * particleCrossSection  * dE)
    neutronIntegral = neutronIntegral + neutronFlux(radius, inclination, E)*neutronCrossSection(radius, inclination, E)*deltaE;   
  }
  // multiplies by the number of targets we have.
  neutronIntegral = neutronIntegral * neutronTargetNumber;
  

  double integral = neutronIntegral;// so I can add other things like protons as I go
  
  return integral;
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
  double Emax = pow(10,7)*kUnits::MeV;
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
  double Emax = pow(10,2)*kUnits::MeV;
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
  
  double background = area * 3.3*pow(10,4)*1.2; // (events/s m^2) * m^2// 120% of the approximate proton rate to account for alphas 
    //integral; // Goes to 0 at 0, and 1 far from the sun
  
  cout << "cosmicBackground has been called and is returning " << background << " events per second.\n";
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


int nuSolPerformanceElipticalMinutes(){
  // These are values that the main program manipulates
  double time, radius, inclination,nuAccumulator,cosmicAccumulator,solarAccumulator=1;
  gRandom -> SetSeed(0);
  double nuThresh = (gRandom->Uniform());
  double cosmicThresh = (gRandom->Uniform());
  double solarThresh = (gRandom->Uniform());
  bool isNeutrino, isCosmic, isSolar, isRadio = false; 
  const double cosBack = cosmicBackground(1,1); // currently flat, so we save time
  cout << "cosmic background is " << cosBack << endl;


  // Control Variables!
  double timeStep = 1000*kUnits::min;// time step for iterating
  string myFile = "solar_distance_C3_60(days).csv"; // String that holds r vs time
  double timeLimit = 1*kUnits::day;//1829*kUnits::day; // constant cap
  double cosmicAcceptance = 0.0001; // 4 nines
  double solarAcceptance = cosmicAcceptance/10;// one nine better than cosmic
  //double timeLimit = number_of_lines(myFile)*kUnits::day; // soft cap

  // These are the histograms we output

  TH1D* neutrinoRadiusHistogram = new TH1D("neutrinoRadiusHistogram","Neutrino Event Radius",150,3,220);// should always be larger than 0, this program should never have to deal with R>~1.1AU~=230 Rsun
  TH1D* sub35NeutrinoRadiusHistogram = new TH1D("sub35NeutrinoRadiusHistogram","Neutrino Event Radius",150,3,38);// should always be larger than 0, this program should never have to deal with R>~1.1AU~=230 Rsun

  TH1D* neutrinoOnlyRadiusHistogram = new TH1D("neutrinoOnlyRadiusHistogram","Neutrino Only Event Radius",150,3,220);// should always be larger than 0, this program should never have to deal with R>~1.1AU~=230 Rsun
  TH1D* sub35NeutrinoOnlyRadiusHistogram = new TH1D("sub35NeutrinoOnlyRadiusHistogram","Neutrino Only Event Radius",150,3,38);// should always be larger than 0, this program should never have to deal with R>~1.1AU~=230 Rsun
  
  TH1D* cosmicBackRadiusHistogram = new TH1D("cosmicBackRadiusHistogram","Cosmic Background Event Radius",150,3,220);// should always be larger than 0, this program should never have to deal with R>~1.1AU~=230 Rsun
  TH1D* sub35CosmicBackRadiusHistogram = new TH1D("sub35CosmicBackRadiusHistogram","Cosmic Background Event Radius",150,3,38);// should always be larger than 0, this program should never have to deal with R>~1.1AU~=230 Rsun

  TH1D* solarBackRadiusHistogram = new TH1D("solarBackRadiusHistogram","Solar Background Event Radius",150,3,220);// should always be larger than 0, this program should never have to deal with R>~1.1AU~=230 Rsun

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
  double closest = 15*kUnits::solarRadii;
  double furthest = 1.1*kUnits::venusSunDistance;
  double a = (closest+furthest)/2;// semi-major axis
  
  //timeLimit = floor ( ( 2*M_PI*sqrt(a*a*a/( kUnits::G * kUnits::mSun ) ) ) / kUnits::day ) * kUnits::day ;//2*M_PI*sqrt( pow(a,3) / ( kUnits::G*kUnits::mSun ) );
  timeLimit = ( 2*M_PI*sqrt(a*a*a/( kUnits::G * kUnits::mSun ) ) );//2*M_PI*sqrt( pow(a,3) / ( kUnits::G*kUnits::mSun ) );
  
  double* quadraticConstants = elliptical(closest,furthest);
  cout << "elliptical  is done.\n\n";
  for (int i = 0;i<100;i++){
    double c = *(quadraticConstants+3*i);
    double b = *(quadraticConstants+3*i+1);
    double a = *(quadraticConstants+3*i+2);
    //cout << "Day " << (i+1) << " to " << (i+2) << " has constants:\n";
    //cout << "c = " << c;
    //cout << "; b = " << b;
    //cout << "; a = " << a << ".";
    //cout << "Ergo the radius on day " << i+1 << " is " << ( (a*(i+1)*(i+1))+(b*(i+1))+c) << ".\n\n";
  }
  //*/


  //TNtuple* nuNtuple = new TNtuple("nuNtuple","something","Solar_Radii:neutrinoCount:cosmicCount:solarCount");
  int sanityInt = 0;
  for(int i = 0; i < 1; i++ ){
      time = 0;
      do{
	sanityInt += 1;
	ofstream sanityCSV;
	string sanityName = "sanityCheck";
	sanityName += std::to_string(i);
	sanityName += ".txt";
	sanityCSV.open (sanityName, std::ios::app);
	cout << "The file name is \"" << sanityName << "\".\n";
	std::cout << "Before flightPath, time is " << (time/kUnits::days) << " d, radius is " << (radius/kUnits::solarRadii) << " RSol, and inclination is " << inclination <<endl;

	// flightPath(time,radius,inclination);
	radius = (positionMinutes( time, quadraticConstants, timeLimit ));
        sanityCSV << "This is step " << sanityInt << ".\n";
	sanityCSV << radius << ',' << (time/kUnits::hr) << endl ;

	if (radius > 250*kUnits::solarRadii){
	  cout << "Oh no! The radius is " << radius/kUnits::solarRadii << " solar radii, which is much further than Earth. I'm setting it to 500 R_Sol \n";
	  radius = 500*kUnits::solarRadii;
	}
	if (radius < 1) {
	  cout << "Well shoot, we're at " << radius << " solar radii, which is inside the sun. (Or negative. That makes even less sense). I'm setting it to 0.5 R_sol.\n";
	  radius = 0.5*kUnits::solarRadii;
	}
    
	time = time + timeStep;
	// std::cout << "After flightPath, time is " << time << " radius is " << radius << ", and inclination is " << inclination <<endl;

	// The GCR is too big for most timescales directly, but small enough that we don't need it to be
	// I can take the cosmic accumulator, and make it the nearest integer. Once I have that, I can 
	// randomly assign events to bins in the timescale range. A simpler approximation would be to just 
	// assign them to bin 0 through n until I'm out of events. If I do them randomly, then I'll have 
	// to look at the probability of having two hit at once, and the veto performance of that. For
	// the time being, I'll probably stick to even bins, and randomly select a bin for the neutrino.
	bool cosmicBack = false;
	if(cosmicBack){   
	  double background = cosmicBackground(radius,1); // checking background
	  cout << "Total backround is " << background << " events per second.\n";
	  double vetoedBack = background * cosmicAcceptance * cosmicAcceptance;
	  cout << "Vetoed backround is " << vetoedBack << " events per second.\n";
	  double a = timeStep * cosmicBackground(radius,1) * cosmicAcceptance * cosmicAcceptance;//cosmicBackground(radius,1);
	  std::cout << "The Cosmic Background of " << (radius/kUnits::solarRadii)  << " solarRadii from the center of the sun is : " << a << endl;
	  cosmicAccumulator = cosmicAccumulator + a;
	  std::cout << "The Cosmic Background Signal accumulation is " << cosmicAccumulator  << ". " << endl;
	  isCosmic = cosmicAccumulator>cosmicThresh;
	  if(isCosmic){
	    cosmicAccumulator = 0;//cosmicAccumulator - 1;
	    cosmicThresh = (gRandom->Uniform());
	    cosmicBackRadiusHistogram -> Fill( (radius / kUnits::solarRadii) );
	    if(radius < 35*kUnits::solarRadii) sub35CosmicBackRadiusHistogram -> Fill( (radius / kUnits::solarRadii) );
	    //nuNtuple->Fill((radius/kUnits::solarRadii),0,1,0);
	  }
	}
  
	bool solarBack = false;
	if(solarBack){  
	  double b = timeStep * solarBackground(radius,1);
	  std::cout << "The Solar Background of " << (radius/kUnits::solarRadii)  << " solarRadii from the center of the sun is : " << b << endl;
	  solarAccumulator = solarAccumulator + b;
	  std::cout << "The Solar Background Signal accumulation is " << solarAccumulator  << ". " << endl;
	  isSolar = solarAccumulator>solarThresh;
	  if(isSolar){
	    solarAccumulator = 0;//solarAccumulator - 1;
	    solarThresh = (gRandom->Uniform());
	    solarBackRadiusHistogram -> Fill( (radius / kUnits::solarRadii) );
	    //nuNtuple->Fill((radius/kUnits::solarRadii),0,0,1);
	  }
	}
  
	bool neutrinos = false;
	if(neutrinos){  
	  double c = timeStep*neutrinoSignal(radius,1);
	  std::cout << "The Neutrino Signal of " << (radius/kUnits::solarRadii)  << " solarRadii from the center of the sun is : " << c << endl;
	  nuAccumulator = nuAccumulator + c;// this creates a data race. I should make a raceless version that is parallel. It may be fast enough to make parallel faster. Each RNG call is ~5ns
	  std::cout << "The Neutrino Signal accumulation is " << nuAccumulator  << ". " << endl;
	  isNeutrino = nuAccumulator>nuThresh;
	  if(isNeutrino){
	    bool neutrinoOnly = isNeutrino && !( isCosmic || isSolar || isRadio);
	    nuAccumulator = 0;//nuAccumulator-1;
	    nuThresh = (gRandom->Uniform());
	    neutrinoRadiusHistogram -> Fill( (radius / kUnits::solarRadii) );
	    if(radius < 35*kUnits::solarRadii) sub35NeutrinoRadiusHistogram -> Fill( (radius / kUnits::solarRadii) );
	    if(neutrinoOnly) {
	      neutrinoOnlyRadiusHistogram -> Fill( (radius / kUnits::solarRadii) );
	      if(radius < 35*kUnits::solarRadii) sub35NeutrinoOnlyRadiusHistogram -> Fill( (radius / kUnits::solarRadii) );
	    }
	    //nuNtuple->Fill((radius/kUnits::solarRadii),1);
	  }
	}
	radiusHistogram -> Fill( (radius / kUnits::solarRadii) );
    
	std::cout<<endl;
      }
      while(time<timeLimit);
    }
  
  // print out the NTuple
  TFile *outfile = TFile::Open("outfile.root","recreate");
  neutrinoRadiusHistogram -> Write();
  sub35NeutrinoRadiusHistogram -> Write();
  neutrinoOnlyRadiusHistogram -> Write();
  sub35NeutrinoOnlyRadiusHistogram -> Write();
  solarBackRadiusHistogram -> Write();
  cosmicBackRadiusHistogram -> Write();
  sub35CosmicBackRadiusHistogram -> Write();
  radiusHistogram -> Write();
  outfile -> Close();

  return 0;
}


