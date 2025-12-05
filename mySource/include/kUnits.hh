/*
  Jonathan Folkerts Feb 2021
  This the the file for importing my system of units to the program. I'll be maintiaing this as I
  try to get nuSolPerformance cleaned up into one program with subfiles.

*/

#ifndef KUNITS_H
#define KUNITS_H

#include <math.h>
#include <complex>

namespace kUnits
{

// Math Units
  const std::complex<double> i(0.0,1.0);// imaginary unit for ease of use


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
  const double solarRadius = 695700*km;

  const double venusSunDistance = 1.075e8*km;

  const double mercurySunDistance = 4.602e7*km;

  const double AU = 149597870700*m;





  const double barn = 1e-28*m*m;
  const double b = barn;
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

  const double month = 30.42*day;
  const double months = month;

  const double year = 365.24219*day;
  const double years = 365.24219*day;
  const double yr = 365.24219*day;



  // Electric Potential Units

  const double volt = 1;
  const double V = volt;

  const double millivolt = V/1000;
  const double mV = V/1000;

  const double microvolt = mV/1000;
  const double uV = mV/1000;

  const double nanovolt = uV/1000;
  const double nV = uV/1000;

  const double kilovolt = V*1000;
  const double kV = V*1000;

  const double megavolt = kV*1000;
  const double MV = kV*1000;




  // Energy Units
  const double MeV = 1;

  const double keV = MeV/1000;
  
  const double eV = keV/1000;
  const double electronVolt = eV;

  const double GeV = 1000*MeV;

  const double Joule = 6.241509e18*eV;
  const double J = Joule;

  // Power Units
  const double Watt = J/s;
  const double Watts = J/s;
  const double W = J/s;
  
  const double kilowatt = 1000*W;
  const double kilowatts = 1000*W;
  const double kW = 1000*W;
  
  const double megawatt = 1000*kW;
  const double megawatts = 1000*kW;
  const double MW = 1000*kW;
  
  const double gigawatt = 1000*MW;
  const double gigawatts = 1000*MW;
  const double GW = 1000*MW;
  
  const double solarLuminosity = 3.828e26*Watt;
  
  

  // Mass Units
  const double kg = s*s*J/m;
  const double kilogram = kg;
  
  const double gram = kg/1000;
  
  const double amu = gram/6.02214076e23;

  const double mSun = 1.9891e30*kg; // might want to swap with kg for 1

  const double me = 548.579909065e-6*amu; // electron mass
  const double electronMass = me;

  // Charge
  const double electronCharge = 1;
  const double elementaryCharge = electronCharge;
  const double e = electronCharge;

  const double coulomb = e/1.602176634e-19;


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
  
  // Temperature Units
  const double Kelvin = 1;
  const double K = Kelvin;
  
  // Physical Constants
  const double c = 299792458*m/s;
  const double speedOfLight = c;

  const double h = 6.62606896e-34*Joule*s;
  //const double h = 4.136e-15*eV*s;
  const double hbar = h/(2*M_PI);

  const double G = 6.67430e-11*m*m*m/kg/s/s; 
  const double newtonGravity = G; 

  const double k_b = 1.380649e-23*kUnits::J/kUnits::K;
  const double boltzmannConstant = k_b;

  const double fermiCouplingConstantNaturalUnits = 1.1663787e-5/GeV/GeV;
  const double fermiCouplingConstant = fermiCouplingConstantNaturalUnits *(hbar*c)*(hbar*c)*(hbar*c);

  const double mol = 6.02214076e23;
  const double N_A = mol;
  const double AvogadroNumber = mol;

  const double vacuumPermittivity = 8.8541878188e-12 * coulomb * coulomb / kg / m / m / m * s * s;
  const double epsilonNought = vacuumPermittivity;
    
  // best fit values from pdglive normla ordering
  const double theta12 = asin(sqrt(0.307));
  const double theta13 = asin(sqrt(0.0220));
  const double theta23 = asin(sqrt(0.546));
  const double deltam21 = 7.53e-5*kUnits::eV*kUnits::eV;
  const double deltam32 = 0.002453*kUnits::eV*kUnits::eV;
  const double deltacp = 1.36;
    
  // best fit values from pdglive normla ordering
  const double theta12IO = asin(sqrt(0.307));
  const double theta13IO = asin(sqrt(0.0238));
  const double theta23IO = asin(sqrt(0.539));
  const double deltam21IO = 7.53e-5*kUnits::eV*kUnits::eV;
  const double deltam32IO = -0.002524*kUnits::eV*kUnits::eV;
  const double deltacpIO = 1.57;


}

#endif // kUnits_H
