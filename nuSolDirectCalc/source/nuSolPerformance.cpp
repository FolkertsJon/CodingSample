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
#include "kUnits.hh"
#include "quadFit.hh"
#include "galliumInteraction.hh"
#include "backgrounds.hh"
#include "interp.hh"
#include "variableParser.hh"
#include "trajectoryToCoords.hh"
using namespace std;


// This is the main program


// too easy to make a data race in this code
//int nThreads = 12;
//double deltaTimeLimit = timeLimit/nThreads;


int main(int argc, char * argv[]){
  // These are values that the main program manipulates
  gRandom -> SetSeed(0);
  bool isNeutrino, isCosmic, isSolar, isRadio = false; 
  string nameOfOutfile = "outfile.root";


  // Control Variables!
  User::Par *variables = new User::Par; // holds control variables details in variableParser.hh
  // parse arguments
  User::argumentInterpreter(argc, argv, variables);
  User::excitedOnly = variables->excitedOnly;
  // Check
  std::cout << "File Mode is: " << variables->fileMode << ".\n";
  std::cout << "The number of loops defaults to 100 and is actually " << variables->nLoops << ".\n";
  std::cout << "Excited state boolean is " << variables->excitedOnly << ".\n";
  std::cout << "The name of the file I'm supposed to open is " << std::string(variables->myFile) << '\n';
  std::cout << "The closest approach is: " << variables->closest/kUnits::solarRadii << " solar radii.\n";
  std::cout << "The furthest approach is: " << variables->furthest/kUnits::solarRadii << " solar radii.\n";
  std::cout << "The number of kilograms I'm simulating is: " << variables->kgGallium/kUnits::kg;
  std::cout << "\nThe gallium atom number I'm simulating is: " << variables->galliumAtomNumber;
  //std::cout << "\nThe time step is: " << variables->timeStep/kUnits::second << "s.\n";
  variables->solarNeutrinoRate = User::neutrinoRate()*variables->galliumAtomNumber*0.399; // update rate
  std::cout << "The rate at earth per " << variables->kgGallium/kUnits::kg << " kg of natural gallium is "
	    << variables->solarNeutrinoRate/kUnits::s << " per second.\n";
  //std::cout << "\nPress return to continue.\n";
  //std::cin.ignore();

  
 
  double a = (variables->closest+variables->furthest)/2;// semi-major axis

  // create outfile name
  nameOfOutfile = "";
  nameOfOutfile += std::to_string(variables->nLoops);
  nameOfOutfile += "Orbits_";
  nameOfOutfile += std::to_string(variables->kgGallium/kUnits::kg);
  std::cout << "kilograms are " << std::to_string(variables->kgGallium/kUnits::kg) << "\n.";
  nameOfOutfile += "kg_VariableTimeStep_";
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
  TNtuple *myTuple = new TNtuple("myTuple","NeutrinoFlightPath","neutrinoRadius:sub35NeutrinoRadius:radius:neutrinoSignal:timeStep");

  double yearFactor;
  
  if(variables->fileMode){
    yearFactor = 1; // no weighting, we want the whole orbit
    double *t = User::trajectoryTot(std::string(variables->myFile));
    variables->timeLimit = kUnits::day*t[User::number_of_lines(std::string(variables->myFile))-1];//kUnits::hours*User::number_of_lines(variables.myFile);
    std::cout << "The time limit is :" <<  variables->timeLimit/kUnits::day << " days\n";
    delete t;
    variables->nSteps = floor(variables->timeLimit/variables->timeStep);
    //variables->nSteps = 10000000; // does this work?
    std::cout << "The number of steps is: " << variables->nSteps << "\n";
  }
  else
    {
      variables->timeLimit = ( 2*M_PI*sqrt(a*a*a/( kUnits::G * kUnits::mSun ) ) );
      // weight so that one orbit becomes 1 year of data
      yearFactor = kUnits::year/variables->timeLimit;
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
    cout << "Or here?";


    for(int j = 0; j < variables->nSteps; j++){
      //std::cout << "The absolute first line is a comment in this for loop.\n";
      // holds the interpolation data
      double tee[5];
      double arr[5];
      
      // find the index of the time step greater than this
      int idx = 0;
      while( j*variables->timeStep/kUnits::day > t[idx]){
	idx++;
      }//end While
      //std::cout << "here. Why? The index is " << idx<< "\n";
      if(idx < 3 ){ // index is smoll
	for(int k = 0; k < 5; k++){
	  tee[k]=t[k];
	  arr[k]=R[k];
	}
      }// end if1
      else if(idx > (i - 3) ){//index is big
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
      output = interp_lagrange( 1, 5, tee, arr, 1, interpPoint);
      fileConstants[j] = output[0];
      //std::cout << "I'm outputting!\n";
      
      
      if(j%10000 ==0){
	std::cout << "On fit value " << j << " of " << variables->nSteps << " and the radius is " << output[0] << ".\n";
      }
      //std::cout << "The absolute last line is a comment in this for loop.\n";
    }//end for
  }//end file mode if
  else{
    std::cout << "The closest approach is: " << variables->closest/kUnits::solarRadii << " solar radii.\n";
    std::cout << "The furthest approach is: " << variables->furthest/kUnits::solarRadii << " solar radii.\n";
    quadraticConstants = User::elliptical(variables->closest,variables->furthest);
  }
  cout << "fitting is done.\n\n";


  for(int i = 0; i < variables->nLoops; i++ ){
    double maxRadiusStep = 1e-6*kUnits::AU; // 10 times smaller than precision of Parker Probe
    double maxTimeStep = 60000*kUnits::sec;
    double time, radius = 0;
    time = 0;
    long long int sanityInt = 0;
    cout << "\n\nThis is iteration " << i << ".\n";
    do{
      //timeLimit = 1;
      sanityInt += 1;
      int coutPeriod;
      coutPeriod = (int)1e6;
      bool doCout = (0 == (sanityInt%coutPeriod) );
      // this is to suppress a NaN issue. This makes the elliptical loops less good, but the
      // explicit file loops are unaffected
      
      //double nuAccumulator,cosmicAccumulator,solarAccumulator=0;
      if(doCout){
	std::cout << "Before flightPath, time is " << (time/kUnits::days) << " d, radius is " << (radius/kUnits::solarRadii) << " RSol, and inclination is " << variables->inclination <<endl;
      }


      // flightPath(time,radius,inclination);
      if(variables->fileMode){
	int thisStep = floor(time/variables->timeStep);
	radius = fileConstants[thisStep]*kUnits::solarRadii;
	//std::cout << "This is step " << thisStep << ".\n";
	//radius = (User::positionHours( time, quadraticConstants, timeLimit ));
      }
      else{
	//std::cout << "Derivative\n";
	variables->timeStep =min(maxRadiusStep/User::derivativeMinutes( time, quadraticConstants, variables->timeLimit ),
				 maxTimeStep);
	//std::cout << "Derivative done\n";
	radius = (User::positionMinutes( time, quadraticConstants, variables->timeLimit ));
      }
      if(std::isnan(radius)){
	radius = 0;// controls a NaN radius
      }
      double solarNeutrinoSignal = variables->solarNeutrinoRate*(1*kUnits::AU/radius)*(1*kUnits::AU/radius);
      if( solarNeutrinoSignal ==  0) std::cout << "Oh god, dear god help. The singal is zero. The signal should never be zero!!!\n\n";
      
      if (radius > 250*kUnits::solarRadii){
	std::cout << "Before flightPath, time is " << (time/kUnits::days) << " d, radius is " << (radius/kUnits::solarRadii) << " RSol, and inclination is " << variables->inclination <<endl;
	cout << "Oh no! The radius is " << radius/kUnits::solarRadii << " solar radii, which is much further than Earth. I'm setting it to 500 R_Sol or breaking.\n";
	radius = 500*kUnits::solarRadii;
	break;
      }
      if (radius < 1*kUnits::solarRadii) {
	std::cout << "Before flightPath, time is " << (time/kUnits::days) << " d, radius is " << (radius/kUnits::solarRadii) << " RSol, and inclination is " << variables->inclination <<endl;
	cout << "Well shoot, we're at " << radius << " solar radii, which is inside the sun. (Or negative. That makes even less sense). I'm setting it to 0.5 R_sol or breaking.\n";
	radius = 0.5*kUnits::solarRadii;
	break;
      }
      
      time = time + variables->timeStep;
      // std::cout << "After flightPath, time is " << time << " radius is " << radius << ", and inclination is " << inclination <<endl;
      
      // The GCR is too big for most timescales directly, but small enough that we don't need it to be
      // I can take the cosmic accumulator, and make it the nearest integer. Once I have that, I can 
      // randomly assign events to bins in the timescale range. A simpler approximation would be to just 
      // assign them to bin 0 through n until I'm out of events. If I do them randomly, then I'll have 
      // to look at the probability of having two hit at once, and the veto performance of that. For
      // the time being, I'll probably stick to even bins, and randomly select a bin for the neutrino.
      
      double c = variables->timeStep*solarNeutrinoSignal;
      if(std::isnan(c)){
	std::cout << "I have somehow made c into NaN. c is just timeStep (" << variables->timeStep 
		  << ") + (" << solarNeutrinoSignal <<").\n"; // solarNeutrinoSignal is somehow NaN
	doCout = true;
	//break;
      }
      double nuAccumulator = c;// this creates a data race. I should make a raceless version that is parallel. It may be fast enough to make parallel faster. Each RNG call is ~5ns

      if(doCout){
	std::cout << "The signal before time and radius stuff is " << solarNeutrinoSignal << std::endl;
	std::cout << "The Neutrino Signal of " << (radius/kUnits::solarRadii)  << " solarRadii from the center of the sun is : " << c << endl;
	std::cout << "The Neutrino Signal accumulation is " << nuAccumulator  << ". " << endl;      }

      {
	
	theHistogramData -> radius = radius / kUnits::solarRadii;
        
	
        
        
	theHistogramData -> neutrinoRadius = radius / kUnits::solarRadii;

	if(radius / kUnits::solarRadii < 35){
	  theHistogramData -> sub35NeutrinoRadius = radius / kUnits::solarRadii;
	}
	  
        

	
	
	//radiusHistogram -> Fill( (radius / kUnits::solarRadii) );
	if(doCout){
	  std::cout << "Filling with: "
		    << theHistogramData -> neutrinoRadius << ','
		    << theHistogramData -> sub35NeutrinoRadius << ','
		    << theHistogramData -> radius << ','
		    << nuAccumulator << ','
		    << variables -> timeStep << endl;
	  std::cout<<endl;

	}
	myTuple ->Fill(theHistogramData -> neutrinoRadius,
	               theHistogramData -> sub35NeutrinoRadius,
	               theHistogramData -> radius,
		       nuAccumulator*yearFactor,
		       variables -> timeStep);
	
	//myTree -> Fill();
      }
      //if (radius < 15.5*kUnits::solarRadii) break;
      /*if(std::isnan(nuAccumulator)||std::isnan(cosmicAccumulator)){
	break;
	}*/
    }
    while(time < variables->timeLimit);
  }

  cout << "\n\nThe total orbital period  was " << variables->timeLimit/kUnits::days << " days.\n\n";
  
  // print out the NTuple
  //myTuple->Draw("neutrinoRadius");
  cout << "I'm writing the histograms now.\n";
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

  TCanvas *c0 = new TCanvas("c0","c0",1920,1080);
  double precision = 1e-6*kUnits::AU/kUnits::solarRadii;
  double binWidth = 1e-6*kUnits::AU/kUnits::solarRadii;
  //double binWidth = 1e-1;
  double xMax = 35 + precision/2;
  double xMin = xMax;
  int nBins = 1;
  while(xMin > 0){
    nBins++;
    xMin -= binWidth;
  }
  
  string title;
  if(variables->fileMode){
    yearFactor = 1;// no weighting
    title = "Number of Neutrinos During ACO Mission Orbit;\\text{Radius(}R_\\odot\\text{)};\\text{Neutrinos (}R_\\odot^{-1} 10y^{-1}\\text{)}";
  }
  else if(variables->furthest/kUnits::solarRadii > 100){
    title = "Number of Neutrinos per Year for Venus to " + std::to_string(variables->closest/kUnits::solarRadii) + " Solar Radii Orbit;\\text{Radius(}R_\\odot\\text{)};\\text{Neutrinos (}R_\\odot^{-1} y^{-1}\\text{)}";
  }
  else{
    title = "Number of Neutrinos per Year for Venus to " + std::to_string(variables->closest/kUnits::solarRadii) + " Solar Radii Orbit;\\text{Radius(}R_\\odot\\text{)};\\text{Neutrinos (}R_\\odot^{-1} y^{-1}\\text{)}";
  }
  TH1D *neutrinoWeighted = new TH1D("neutrinoWeighted",title.c_str(),nBins,xMin,xMax);

  Long64_t nEntries = myTuple -> GetEntries();
  
  float radius, weight;
  myTuple -> SetBranchAddress("radius",&radius);
  myTuple -> SetBranchAddress("neutrinoSignal",&weight);
  
  for(long long int i = 1; i <= nEntries; i++){
    myTuple -> GetEntry(i);
    if(i % 100000 == 0){
      std::cout << "Event " << i << " of " << nEntries << "\n";
      std::cout << "Radius = " << radius << " with weight " << weight << "\n";
    }
    if(radius > 1 && radius < 35){
      //neutrinoWeighted -> Fill(radius, weight*yearFactor);
      neutrinoWeighted -> Fill(radius, weight);
    }
  }

  c0 -> SetLogy();
  neutrinoWeighted -> Draw();

  c0 -> Print("radiusHist.png");



  
  myTuple -> Print();
  myTuple -> Write();
  neutrinoWeighted->Write();
  //myTree -> Print();
  //myTree -> Write();
  outfile -> Close();
  
  
  
  return 0;
}


