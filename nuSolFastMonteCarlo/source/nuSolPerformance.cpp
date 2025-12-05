/*
  This is a first pass at a program designed to perform a monte carlo
  simulation of the nuSol space probe's performance. The first 
  version is mean to be minimally functional, and will include wildly 
  inaccurate models.

  We generate a random number in [0,2] (for an averate of 1). The 
  simulation will begin stepping through the flight path until the 
  probability of seeing a neutrino event exceeds this random number. 
  We will take this event and give it energy and direction. Then we 
  will smear these values to represent our uncertainties. We apply 
  an effeciency module and then write this event to a ROOT tree. 
  This continues until the detector reaches the end of the flight path.


  Version history (Should ahve been doing this more)
  - Jan 19: Exits
  Jan 19 2021 : Removed /steradian from the neutrino rates. (units in graph 
                they were taken from did not have them. Oops)
  Feb 2       : Removing code into sub-files: kUnits,
  Feb 4       : sub-files: quadFit (does the elliptical or file quadratic 
                fitting)
  May 11      : Program now can interpolate non-linear positions; this is 
                held in a large double array, which might be problmatic.
		The loops don't seem to be working properly. The neutrino
		counts don't get looped more than once by appearance
  May 8 2023  : some time earlier, I added updated fluxes from the Bahcall;
                also I'm updating the histograms to weight their entries by
		1/number of flights; jk, ntuples can't be weighted

  
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
#include <TTree.h>
#include <TNtuple.h>
#include <TH1D.h>
#include <TRandom.h>
#include "TParameter.h"

#include "kUnits.hh"
#include "quadFit.hh"
#include "galliumInteraction.hh"
#include "backgrounds.hh"
#include "interp.hh"
#include "variableParser.hh"
#include "trajectoryToCoords.hh"


// This is the main program


// too easy to make a data race in this code
//int nThreads = 12;
//double deltaTimeLimit = timeLimit/nThreads;


int main(int argc, char * argv[]){
  if (argc % 2 == 0){
    User::printHelp();
    return 0;
  }

  
  // These are values that the main program manipulates
  gRandom -> SetSeed(0);
  bool isNeutrino, isCosmic, isSolar, isRadio = false; 
  std::string nameOfOutfile = "outfile.root";// default name for output file


  // Control Variables!
  User::Par *variables = new User::Par; // holds control variables details in variableParser.hh
  User::doOscillations = false;
  // parse arguments
  
  User::argumentInterpreter(argc, argv, variables);
  if(User::printedHelp) return 0;
  User::excitedOnly = variables->excitedOnly; // false unless set in macro

  
  std::cout << "excitedOnly defaults to false and is: " << User::excitedOnly << "\n";
  std::cout << "fileMode is set to: " << variables->fileMode << "\n";
  // Check
  std::cout << "The number of loops defaults to 100 and is actually " << variables->nLoops << ".\n";
  std::cout << "The time limit defaults to 0  and is actually " << variables->timeLimit/kUnits::year << " years.\n";
  std::cout << "The name of the file I'm supposed to open is " << std::string(variables->myFile) << '\n';
  std::cout << "The closest approach is: "<< variables->closest/kUnits::solarRadii << " solar radii.\n";
  std::cout << "The furthest approach is: " << variables->furthest/kUnits::solarRadii << " solar radii.\n";
  std::cout << "The number of kilograms I'm simulating is: " << variables->kgGallium/kUnits::kg;
  std::cout << "\nThe gallium atom number I'm simulating is: " << variables->galliumAtomNumber;
  std::cout << "\nThe time step is: " << variables->timeStep/kUnits::second << "s.\n";

  std::cout << "Old per atom rate: ";
  variables->solarNeutrinoRate = User::neutrinoRate(0,20)*variables->galliumAtomNumber*0.399; // update rate
  std::cout << "* " << variables->galliumAtomNumber << " atoms\n"
	    << "* 0.399 \ngives:" << variables->galliumAtomNumber*0.399 << " atoms\n";
    std::cout << "\nThe updated neutrino rate for the current mass of natural gallium is: " << variables->solarNeutrinoRate*kUnits::s << " v/s.\n";
  

  // to pause here if need be for checks
  //std::cout << "\nPress return to continue.\n";
  //std::cin.ignore();

  
  
  // Prepare file names and orbital parameters
  double a = (variables->closest+variables->furthest)/2;// semi-major axis

  // create outfile name
  nameOfOutfile = "";
  if(variables->timeLimit > 0){
    nameOfOutfile += std::to_string(variables->timeLimit/kUnits::yr);
    nameOfOutfile += "years_";    
  }
  else{
    nameOfOutfile += std::to_string(variables->nLoops);
    nameOfOutfile += "Orbits_";
  }
  nameOfOutfile += std::to_string(variables->kgGallium/kUnits::kg);
  std::cout << "kilograms are " << std::to_string(variables->kgGallium/kUnits::kg) << ".\n";
  nameOfOutfile += "kg_";
  nameOfOutfile += std::to_string(variables->timeStep/kUnits::sec);
  nameOfOutfile += "secondTimeStep_";
  if(variables->fileMode){
    nameOfOutfile += "FileOrbit";
  }
  else{
    nameOfOutfile += "EllipticalOrbit_Closest=";
    nameOfOutfile += std::to_string(variables->closest/kUnits::solarRadii);
    nameOfOutfile += "_Furthest=";
    nameOfOutfile += std::to_string(variables->furthest/kUnits::solarRadii);
  }
  nameOfOutfile += ".root";
  
  // This is the outfile
  const char* nameHolder = nameOfOutfile.c_str();
  TFile *outfile = TFile::Open(nameHolder,"recreate");
  


  // These are the histograms we output
  User::HistoDat* theHistogramData = new User::HistoDat; // structure to hold data
  //User::HistoDat theHistogramData;
  TNtuple *myTuple = new TNtuple("myTuple","NeutrinoFlightPath","neutrinoRadius:sub35NeutrinoRadius:cosmicBackRadius:sub35CosmicBackRadius:solarBackRadius:sub35SolarBackRadius");
  TNtuple *radiusTuple = new TNtuple("radiusTuple","NeutrinoFlightPath","radius");
  Double_t fractionalNuSeen = 0;
  Double_t fractionalNuThere = 0;

  //
  
  std::vector<std::vector<double>> RadiusAndTime;
  TGraph *RadiusVsTime;

  
  double orbitTimeLimit = 0;
  if(variables->fileMode){
    double *t = User::trajectoryTot(std::string(variables->myFile));
    orbitTimeLimit = kUnits::day*t[User::number_of_lines(std::string(variables->myFile))-1];
    std::cout << "The orbit time is: " <<  orbitTimeLimit/kUnits::day << " days\n";
    delete t;
    variables->nSteps = floor(orbitTimeLimit/variables->timeStep);
    std::cout << "The number of steps is: " << variables->nSteps << "\n";
  }
  else{
    orbitTimeLimit = ( 2*M_PI*sqrt(a*a*a/( kUnits::G * kUnits::mSun ) ) );
    std::cout << "The orbit time is: " <<  orbitTimeLimit/kUnits::day << " days\n";
  }
  
  // Make Constants the right size
  double* fileConstants = new double[variables->nSteps];
  double* quadraticConstants;
  if(variables->fileMode){
    // Get array size
    int i = User::number_of_lines(std::string(variables->myFile));
    
    //Get Data
    double *R = new double[i];
    double *t = new double[i];
    R = User::trajectoryToR(std::string(variables->myFile));
    t = User::trajectoryTot(std::string(variables->myFile));
    

    for(int j = 0; j < variables->nSteps; j++){
      // holds the interpolation data
      double tee[5]; // time 
      double arr[5]; // radius
      
      // find the index of the time step greater than this
      int idx = 0;
      while( j*variables->timeStep/kUnits::day > t[idx]){
	idx++;
      }//end While
      if(idx < 3 ){ // index is smol
	for(int k = 0; k < 5; k++){
	  tee[k]=t[k];
	  arr[k]=R[k];
	}
      }// end if1
      else if(idx > (i - 3) ){// index is big
	int kay = 0;
	for(int k = variables->nSteps-5; k < variables->nSteps; k++){
	  tee[kay]=t[k];
	  arr[kay]=R[k];
	  kay++;
	}
      }//end if2
      else{
	int kay = 0;
	for(int k = idx-2; k < idx+3; k++){
	  //std::cout << "here. Why? k = "<< k << ". The value of t is " << t[k] << ".\n";
	  tee[kay]=t[k];
	  arr[kay]=R[k];
	  kay++;
	}
      }// end else


      
      //time to interpolate
      double theValue = j*variables->timeStep/kUnits::day;
      double interpPoint[1] = {theValue};
      double *output = new double[1];
      output = interp_lagrange( 1, 5, tee, arr, 1, interpPoint); // fancy interpolation from library
      fileConstants[j] = output[0];
      //std::cout << "I'm outputting!\n";
      
      int modNum = int(pow(10,floor(log10(variables->nSteps)) -  2)); // mod number is between 1/100th and 1/1000th of the total number of points
      
      if(j%modNum == 0){
	std::cout << "On fit value " << j << " of " << variables->nSteps << " and the radius is " << output[0] << ".\n";
      }
    }//end for
  }//end file mode if
  else{
    std::cout << "The closest approach is: " << variables->closest/kUnits::solarRadii << " solar radii.\n";
    std::cout << "The furthest approach is: " << variables->furthest/kUnits::solarRadii << " solar radii.\n";
    quadraticConstants = User::elliptical(variables->closest,variables->furthest);
  }
  std::cout << "fitting is done.\n\n";
  

  double nuAccumulator = 0, cosmicAccumulator = 0, solarAccumulator = 0;
  double nuThresh = (gRandom->Exp(1));
  double cosmicThresh = (gRandom->Exp(1));
  double solarThresh = (gRandom->Exp(1));

  
  double time;
  if(variables->timeLimit > 0  && !variables->fileMode){
    // start at a random time in the elliptical orbit to not favor any part of
    // the orbit due to starting at t = 0 (furthest approach)
    time = (gRandom->Uniform(orbitTimeLimit));
  }
  else{
    // time = 0 is fine for n-many complete loops
    time = 0;
  }
  double globalTime = 0;
  size_t i = 0;

  size_t neutrinoCount = 0;
  size_t neutrinoCount35 = 0;

  // hasn't been a problem because the timed runs have never hit 100 loops; this needs fixing
  while(i < variables->nLoops){ // while 1
    if((variables->timeLimit <= globalTime) && variables->timeLimit > 0) {
      std::cout << "Breaking while loop due to timeLimit reached at t = " << variables->timeLimit/kUnits::day << " d\n";
      break;
    }
    double radius = 0;
    long long int sanityInt = 0;
    if( i > 0){
      time = 0;
    }
    std::cout << "\n\nThis is iteration " << i << ".\n";
    do{ // do for while 2
      // break at time limit; no message because while 1 will take care of that
      if((variables->timeLimit <= globalTime) && variables->timeLimit > 0) break;
      //timeLimit = 1;
      int coutPeriod;
      coutPeriod = floor((orbitTimeLimit)/(4*(variables->timeStep)));
      bool doCout = (0 == (sanityInt%coutPeriod) );
      sanityInt += 1;
      // this is to suppress a NaN issue. This makes the elliptical loops less good, but the
      // explicit file loops are unaffected
      
      //double nuAccumulator,cosmicAccumulator,solarAccumulator=0;
      

      /* // CSV is slow
	 ofstream sanityCSV;
	 std::string sanityName = "sanityCheck";
	 sanityName += std::to_string(i);
	 sanityName += ".txt";
	 sanityCSV.open (sanityName, std::ios::app);
	 std::cout << "The file name is \"" << sanityName << "\".\n";
      */
      if(doCout){
	std::cout << "Before flightPath, time is " << (time/kUnits::days) << " d, global time is " << globalTime/kUnits::day <<"d , radius is " << (radius/kUnits::solarRadii) << " RSol, and inclination is " << variables->inclination <<std::endl;
      }


      // flightPath(time,radius,inclination);
      if(variables->fileMode){
	int thisStep = floor(time/variables->timeStep);
	radius = fileConstants[thisStep]*kUnits::solarRadii;
	//std::cout << "This is step " << thisStep << ".\n";
	//radius = (User::positionHours( time, quadraticConstants, timeLimit ));
      }
      else{
	//radius = (User::positionMinutes( time, quadraticConstants, orbitTimeLimit ));
	radius = User::ellipticalRadius(variables->closest, variables->furthest, time); // calculates directly the radius at that time
      }
      if(std::isnan(radius)){
	radius = 0;// controls a NaN radius
      }
      double solarNeutrinoSignal = variables->solarNeutrinoRate*(1*kUnits::AU/radius)*(1*kUnits::AU/radius);
      // with linear increase model
      //double solarNeutrinoSignal = variables->solarNeutrinoRate*User::linOscFactor(radius)*(1*kUnits::AU/radius)*(1*kUnits::AU/radius);
      /*
	std::cout << "The file name is \"" << sanityName << "\".\n";
	std::cout << "Before flightPath, time is " << (time/kUnits::days) << " d, radius is " << (radius/kUnits::solarRadii) << " RSol, and inclination is " << inclination <<std::endl;
	sanityCSV << "This is step " << sanityInt << ".\n";
	sanityCSV << radius << ',' << (time/kUnits::hr) << std::endl ;
      */ // writing to csv is slow
      
      if (radius > 250*kUnits::solarRadii){
	std::cout << "Before flightPath, time is " << (time/kUnits::days) << " d, radius is " << (radius/kUnits::solarRadii) << " RSol, and inclination is " << variables->inclination <<std::endl;
	std::cout << "Oh no! The radius is " << radius/kUnits::solarRadii << " solar radii, which is much further than Earth. I'm breaking.\n";
	radius = 500*kUnits::solarRadii;
	break;
      }
      if (radius < 1*kUnits::solarRadii) {
	std::cout << "Before flightPath, time is " << (time/kUnits::days) << " d, radius is " << (radius/kUnits::solarRadii) << " RSol, and inclination is " << variables->inclination <<std::endl;
	std::cout << "Well shoot, we're at " << radius/kUnits::solarRadii << " solar radii, which is inside the sun. (Or negative. That makes even less sense). I'm breaking.\n";
	radius = 0.5*kUnits::solarRadii;
	break;
      }
      
      time += variables->timeStep;
      globalTime +=  variables->timeStep;
      // std::cout << "After flightPath, time is " << time << " radius is " << radius << ", and inclination is " << inclination <<std::endl;
      
      // The GCR is too big for most timescales directly, but small enough that we don't need it to be
      // I can take the cosmic accumulator, and make it the nearest integer. Once I have that, I can 
      // randomly assign events to bins in the timescale range. A simpler approximation would be to just 
      // assign them to bin 0 through n until I'm out of events. If I do them randomly, then I'll have 
      // to look at the probability of having two hit at once, and the veto performance of that. For
      // the time being, I'll probably stick to even bins, and randomly select a bin for the neutrino.
      bool cosmicBack = true;
      if(cosmicBack){   
	double background = User::cosmicBackground(radius,variables->inclination); // checking background
	double vetoedBack = background * variables->cosmicAcceptance * variables->cosmicAcceptance;
	double a = variables->timeStep * /*cosmicBackground(radius,1)*/ variables->cosBack * variables->cosmicAcceptance * variables->cosmicAcceptance;//cosmicBackground(radius,1);
	cosmicAccumulator = cosmicAccumulator + a;

	if(doCout){
	  std::cout << "Total cosmic backround is " << background << " events per second.\n";
	  std::cout << "Vetoed cosmic backround is " << vetoedBack << " events per second.\n";
	  std::cout << "The Cosmic Background of " << (radius/kUnits::solarRadii)  << " solarRadii from the center of the sun is : " << a << std::endl;
	  std::cout << "The Cosmic Background Signal accumulation is " << cosmicAccumulator  << ". " << std::endl;
	}

	isCosmic = cosmicAccumulator>cosmicThresh;
      }
      
      bool solarBack = true;
      if(solarBack){  
	double background = User::solarBackground( radius, variables->inclination, 2*variables->detectorRadius ); // checking background
	double vetoedBack = background * variables->solarAcceptance * variables->solarAcceptance;
	double b = variables->timeStep * User::solarBackground( radius, variables->inclination, 2*variables->detectorRadius )* variables->solarAcceptance * variables->solarAcceptance;
	solarAccumulator = solarAccumulator + b;

	if(doCout){
	  std::cout << "Total solar backround is " << background << " events per second.\n";
	  std::cout << "Vetoed solar backround is " << vetoedBack << " events per second.\n";
	  std::cout << "The Solar Background of " << (radius/kUnits::solarRadii)  << " solarRadii from the center of the sun is : " << b << std::endl;
	  std::cout << "The Solar Background Signal accumulation is " << solarAccumulator  << ". " << std::endl;
	}

	isSolar = solarAccumulator>solarThresh;
      }
      
      bool neutrinos = true;
      if(neutrinos){  
	double c = variables->timeStep*solarNeutrinoSignal;
	if(std::isnan(c)){
	  std::cout << "I have somehow made c into NaN. c is just timeStep (" << variables->timeStep 
		    << ") * (" << solarNeutrinoSignal <<").\n"; // solarNeutrinoSignal is somehow NaN
	  doCout = true;
	  //break;
	}
	// 4 pi R^2 * flux * AU^2/R^2 = 4 pi AU^2 * flux
	// could be screwed up if there are two neutrinos at one location, but this is vanishingly unlikely
	fractionalNuThere += 4*M_PI*(kUnits::AU * kUnits::AU)*( User::totalFlux() )*(variables->timeStep);
	fractionalNuSeen += c; // just counts the fractions of a neutrino seen, no randomness
	nuAccumulator = nuAccumulator + c;// this creates a data race. I should make a raceless version that is parallel. It may be fast enough to make parallel faster. Each RNG call is ~5ns

	if(doCout){
	  std::cout << "The signal before time and radius stuff is " << solarNeutrinoSignal << std::endl;
	  std::cout << "The Neutrino Signal of " << (radius/kUnits::solarRadii)  << " solarRadii from the center of the sun is : " << c << std::endl;
	  std::cout << "The Neutrino Signal accumulation is " << nuAccumulator  << ". " << std::endl;
	  std::cout<<std::endl;
	}

	isNeutrino = nuAccumulator>nuThresh;
      }
      
      // This is the filling loop
      bool needRad = ( (i == 1) && !(variables->fileMode)) || ( variables->fileMode && (i == 0) );// only need to see radius once; we skip the first partial loop for elliptical
      
      while(isCosmic || isSolar || isNeutrino || needRad){
	if(needRad){
	  theHistogramData -> radius = radius / kUnits::solarRadii;
	  needRad = false;
	}
	else{
	  theHistogramData -> radius = 0;
	}
	
	if(isCosmic){
	  // std::cout << "1/1 There is a cosmic event and the radius I'm writing is: " << radius / kUnits::solarRadii << '\n';
	  theHistogramData -> cosmicBackRadius = radius / kUnits::solarRadii;
	  if(radius / kUnits::solarRadii < 35){
	    theHistogramData -> sub35CosmicBackRadius = radius / kUnits::solarRadii;
	  }
	  cosmicAccumulator = cosmicAccumulator - cosmicThresh;
	  cosmicThresh = (gRandom->Exp(1));
	  isCosmic = (cosmicAccumulator > cosmicThresh);
	  
	}
	else{
	  theHistogramData -> cosmicBackRadius = 0;
	  theHistogramData -> sub35CosmicBackRadius = 0;
	}
      

	if(isSolar){
	  theHistogramData -> solarBackRadius = radius / kUnits::solarRadii;
	  if(radius / kUnits::solarRadii < 35){
	    theHistogramData -> sub35SolarBackRadius = radius / kUnits::solarRadii;
	  }
	  solarAccumulator = solarAccumulator - solarThresh;
	  solarThresh = (gRandom->Exp(1));
	  isSolar = (solarAccumulator > solarThresh);
	  
	}
	else{
	  theHistogramData -> solarBackRadius = 0;
	  theHistogramData -> sub35SolarBackRadius = 0;
	}
      

	if(isNeutrino){
	  neutrinoCount++;
	  theHistogramData -> neutrinoRadius = radius / kUnits::solarRadii;
	  if(radius / kUnits::solarRadii < 35){
	    theHistogramData -> sub35NeutrinoRadius = radius / kUnits::solarRadii;
	  neutrinoCount35++;
	  }
	  nuAccumulator = nuAccumulator - nuThresh;
	  nuThresh = (gRandom->Exp(1));
	  isNeutrino = (nuAccumulator > nuThresh);
	  
	}
	else{
	  theHistogramData -> neutrinoRadius = 0;
	  theHistogramData -> sub35NeutrinoRadius = 0;
	}



	//radiusHistogram -> Fill( (radius / kUnits::solarRadii) );
	// if any accumulators are nonzero
	if(theHistogramData -> neutrinoRadius ||
	   theHistogramData -> sub35NeutrinoRadius ||
	   theHistogramData -> cosmicBackRadius ||
	   theHistogramData -> sub35CosmicBackRadius ||
	   theHistogramData -> solarBackRadius ||
	   theHistogramData -> sub35SolarBackRadius
	   ){
	myTuple ->Fill(theHistogramData -> neutrinoRadius,
	               theHistogramData -> sub35NeutrinoRadius,
	               theHistogramData -> cosmicBackRadius,
		       theHistogramData -> sub35CosmicBackRadius,
	               theHistogramData -> solarBackRadius,
                       theHistogramData -> sub35SolarBackRadius);
	}
	if(theHistogramData -> radius) radiusTuple ->Fill(theHistogramData -> radius);
	
	//myTree -> Fill();
      }
      //if (radius < 15.5*kUnits::solarRadii) break;
      /*if(std::isnan(nuAccumulator)||std::isnan(cosmicAccumulator)){
	break;
	}*/
    } // end do
    while(time < orbitTimeLimit);// while 2
    std::cout << "Finished iteration " << i << ". Beginning " << i+1 <<".\n";
    i++;
  } // end while 1

  std::cout << "\n\nThe total orbital period  was " << orbitTimeLimit/kUnits::days << " days. I saw " << neutrinoCount<< " neutrinos, " << neutrinoCount35 << " of which were inside 35 solar radii.\n\n";

  std::ofstream file("pois.csv", std::ios_base::app);

  file << neutrinoCount35 << '\n';

  file.close();

  
  
  
  // print out the NTuple
  //myTuple->Draw("neutrinoRadius","neutrinoRadius>0");
  std::cout << "I'm writing the histograms now.\n";
  /*neutrinoRadiusHistogram -> Write();
    sub35NeutrinoRadiusHistogram -> Write();
    neutrinoOnlyRadiusHistogram -> Write();
    sub35NeutrinoOnlyRadiusHistogram -> Write();
    solarBackRadiusHistogram -> Write();
    sub35SolarBackRadiusHistogram -> Write();
    cosmicBackRadiusHistogram -> Write();
    sub35CosmicBackRadiusHistogram -> Write();
    radiusHistogram -> Write();
  */
  
  TParameter<Double_t> *fracNuSeen = new TParameter<Double_t>("fractionalNuSeen", fractionalNuSeen);
  TParameter<Double_t> *fracNuThere = new TParameter<Double_t>("fractionalNuThere", fractionalNuThere);
  myTuple -> Print();
  myTuple -> Write();
  radiusTuple -> Write();
  fracNuSeen -> Write();
  fracNuThere -> Write();
  //myTree -> Print();
  //myTree -> Write();
  outfile -> Close();
  std::cout << "\n\n" <<std::flush;
  double neutrinotRate = User::neutrinoRate(0,20);
  std::cout << "Linear expected fluence was " << kUnits::AU*kUnits::AU/kUnits::solarRadius/kUnits::solarRadius*kUnits::day*neutrinotRate << "\n";
  std::cout << "*4*pi*sigma =  " << 4*M_PI*kUnits::AU*kUnits::AU/kUnits::solarRadius/kUnits::solarRadius*kUnits::day*neutrinotRate*variables->galliumAtomNumber*0.399 << "\n";
  
  return 0;
}


