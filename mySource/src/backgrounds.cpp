
// This subprogram will contian the background models for the solar, 
// cosmic, and radiological background processes. These might become
// separate files later, depending on how big this grows. 
// JF - 2/2021



#include "backgrounds.hh"
namespace User
{

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

  double neutronTargetNumber(){
    return 1e27; // 10 times gallium target number for now
      } 
  double protonTargetNumber() {
    return 1e28; // ~assuming the target is 400 pounds of tungsten
  }
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

    //return flux;
    return earthValue *1.5; // assume that up to 1.5x magnification happens, since helios data doesn't show increasing as 1/r^2
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


}// end User namespace












