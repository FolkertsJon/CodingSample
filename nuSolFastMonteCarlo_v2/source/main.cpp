/*
  This is a redesign of the program to do a fast monte carlo 
  of neutrino counts during flights of the nuSOL spacecraft
  I am rewriting this from the ground up to try and get rid 
  of a bug relating to incorrect neutrino counts for 
  non-circular orbits.
  
 */

#include<iostream>
#include<fstream>
#include<cmath>

#include<omp.h>

#include <TFile.h>
#include <TNtupleD.h>
#include <TRandom.h>
#include "TParameter.h"

#include "kUnits.hh"
#include "galliumInteraction.hh"
#include "variableParser.hh"
#include "quadFit.hh"
#include "trajectoryToCoords.hh"



int main(int argc, char * argv[]){
  if (argc % 2 == 0){
    User::printHelp();
    return 0;
  }

  if(User::printedHelp) return 0;
  
  gRandom -> SetSeed(0); // sets seed with current time (0) or a specific seed (nonzero)
  std::string nameOfOutfile = "outfile.root";// default name for output file
  omp_set_num_threads(1);
  
  User::Par *variables = new User::Par; // holds control variables details in variableParser.hh
  
  User::argumentInterpreter(argc, argv, variables);
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


  std::cout << "\nThe expected number of neutrinos for this detector moving 1 solar radius per day from 1 to 100 radii is: "
	    << variables->solarNeutrinoRate*kUnits::AU*kUnits::AU*kUnits::day/kUnits::solarRadius/kUnits::solarRadius*99/100 << "\n";
  

    
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

  
  TFile *outfile = TFile::Open(nameOfOutfile.c_str(),"recreate");

  
  // These are the tuples and doubles we output
  //User::HistoDat theHistogramData;
  TNtupleD *myTuple = new TNtupleD("myTuple","NeutrinoFlightPath","neutrinoRadius:sub35NeutrinoRadius:cosmicBackRadius:sub35CosmicBackRadius:solarBackRadius:sub35SolarBackRadius");
  TNtupleD *radiusTuple = new TNtupleD("radiusTuple","NeutrinoFlightPath","radius");
  Double_t fractionalNuSeen = 0;
  Double_t fractionalNuThere = 0;
  
  std::vector<std::pair<double,double>> orbitPath;
  if(variables->fileMode){
    std::vector<std::pair<double,double>> orbitData = User::radiusTimePairsFromFile(variables->myFile);
    
    double time = 0;
    while(time < orbitData.back().second){
      orbitPath.push_back(std::make_pair(User::interpPosition(orbitData,time),time));
      time += variables->timeStep;
    }
  }
  else{
    std::vector<std::vector<double>> temp = User::ellipticalRadiusTime(variables->closest, variables->furthest,variables->timeStep);
    for(std::vector<double> v : temp){
      orbitPath.push_back(std::make_pair(v[0],v[1]));
    }
  }

  
  // set the number of loops to run if there is a time limit
  std::cout << "Preparing Vector to hold the orbit.\n\n";
  if(variables->timeLimit > 0){
    std::cout << "There is a time limit, so I'm setting the number of loops to ceil( timeLimit/orbitTime)\n\n";
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
    variables->nLoops = ceil(variables->timeLimit/orbitTimeLimit);
    
  }
  std::cout << "Single loop vector is prepared. Orbits begin.\n\n";		      
  

  
  size_t loopCount = 0;
  size_t time = 0;
  double nuAccumulator = 0;
  double nuThresh = (gRandom->Exp(1)); // random on exponential dist

  
  // loop over all but the last loop
  // can be parallalized
#pragma omp parallel for
  for(size_t i = 0; i < (variables->nLoops-1); i++){
    int threadNum = omp_get_thread_num();  // Get the current thread number
    std::cout << "Starting loop " << i << " on thread " << threadNum << ".\n";
    // each thread has its own data
    double fracSeenThisLoop = 0;
    double fracThereThisLoop = 0;
    
    for(std::pair<double,double> radTime : orbitPath){
      User::HistoDat* theHistogramData = new User::HistoDat; // structure to hold data
      bool isNeutrino = false;
      double thisRadius = radTime.first;

      // set the radius to nonzero if on first loop
      if(i == 0){
	theHistogramData -> radius = thisRadius / kUnits::solarRadii;
      }

      // radius independent
      double nuAccumulating = variables->solarNeutrinoRate;
      nuAccumulating *= variables->timeStep;

      // accounts for radius
      nuAccumulating *= kUnits::AU*kUnits::AU/(thisRadius*thisRadius);

      fracSeenThisLoop += nuAccumulating;
      fracThereThisLoop += 4*M_PI*(kUnits::AU * kUnits::AU)*( User::totalFlux() )*(variables->timeStep);
	
#pragma omp critical
      {
	// update global time
	time += variables->timeStep;
	
	nuAccumulator += nuAccumulating;
	// if there's a neutrino to do
	if(nuAccumulator > nuThresh){
	  isNeutrino = true;
	  nuAccumulator -= nuThresh;
	  nuThresh = (gRandom->Exp(1));
	}
      }// end critial

      if(isNeutrino){
	// each thread has its own data
	theHistogramData -> neutrinoRadius = thisRadius / kUnits::solarRadii;
	if(thisRadius / kUnits::solarRadii < 35){
	  theHistogramData -> sub35NeutrinoRadius = thisRadius / kUnits::solarRadii;
	}
      }



      // if any radius values are nonzero
      if(theHistogramData -> neutrinoRadius ||
	 theHistogramData -> sub35NeutrinoRadius ||
	 theHistogramData -> cosmicBackRadius ||
	 theHistogramData -> sub35CosmicBackRadius ||
	 theHistogramData -> solarBackRadius ||
	 theHistogramData -> sub35SolarBackRadius
	 ){
#pragma omp critical
	{
	  // ROOT is not thread safe
	  myTuple ->Fill(theHistogramData -> neutrinoRadius,
			 theHistogramData -> sub35NeutrinoRadius,
			 theHistogramData -> cosmicBackRadius,
			 theHistogramData -> sub35CosmicBackRadius,
			 theHistogramData -> solarBackRadius,
			 theHistogramData -> sub35SolarBackRadius);
	}
      }
      if(theHistogramData -> radius) {
#pragma omp critical
	{
	  Double_t temp[1] = {theHistogramData -> radius};
	  // ROOT is not thread safe
	  radiusTuple ->Fill(temp);
	}
      }
      delete theHistogramData;
    }// end pair for

#pragma omp critical
    {
      fractionalNuSeen += fracSeenThisLoop;
      fractionalNuThere += fracThereThisLoop;
    }
    
  }// end parallel for




  
  // implement last loop here
  // separated because last loop can be not thread safe
  // omp commands are not needed, but easier to leave in
  {
    std::cout << "Starting final loop on one thread.\n";
    double fracSeenThisLoop = 0;
    double fracThereThisLoop = 0;

    // to start at a random point in the orbit
    int rand = gRandom->Integer(orbitPath.size()); // random between 0 and size-1
    
    //std::ofstream outFile("radTimeData.csv");
    for(size_t i = 0; i < orbitPath.size(); i++){
      if(time - variables->timeLimit > 0  && variables->timeLimit > 0) break;
      std::pair<double,double> radTime = orbitPath[(i+rand)%orbitPath.size()];
      //outFile << radTime.first/kUnits::solarRadius << "," << radTime.second/kUnits::sec << "\n";
      User::HistoDat* theHistogramData = new User::HistoDat; // structure to hold data
      bool isNeutrino = false;
      double thisRadius = radTime.first;

      // set the radius to nonzero if this is only loop
      if(variables->nLoops == 1){
	theHistogramData -> radius = thisRadius / kUnits::solarRadii;
      }

      // radius independent
      double nuAccumulating = variables->solarNeutrinoRate;
      nuAccumulating *= variables->timeStep;

      // accounts for radius
      nuAccumulating *= kUnits::AU*kUnits::AU/(thisRadius*thisRadius);

      fracSeenThisLoop += nuAccumulating;
      fracThereThisLoop += 4*M_PI*(kUnits::AU * kUnits::AU)*( User::totalFlux() )*(variables->timeStep);
	
#pragma omp critical
      {
	// update global time
	time += variables->timeStep;
	
	nuAccumulator += nuAccumulating;
	// if there's a neutrino to do
	if(nuAccumulator > nuThresh){
	  isNeutrino = true;
	  nuAccumulator -= nuThresh;
	  nuThresh = (gRandom->Exp(1));
	}
      }// end critial

      if(isNeutrino){
	// each thread has its own data
	theHistogramData -> neutrinoRadius = thisRadius / kUnits::solarRadii;
	if(thisRadius / kUnits::solarRadii < 35){
	  theHistogramData -> sub35NeutrinoRadius = thisRadius / kUnits::solarRadii;
	}
      }



      // if any radius values are nonzero
      if(theHistogramData -> neutrinoRadius ||
	 theHistogramData -> sub35NeutrinoRadius ||
	 theHistogramData -> cosmicBackRadius ||
	 theHistogramData -> sub35CosmicBackRadius ||
	 theHistogramData -> solarBackRadius ||
	 theHistogramData -> sub35SolarBackRadius
	 ){
#pragma omp critical
	{
	  // ROOT is not thread safe
	  myTuple ->Fill(theHistogramData -> neutrinoRadius,
			 theHistogramData -> sub35NeutrinoRadius,
			 theHistogramData -> cosmicBackRadius,
			 theHistogramData -> sub35CosmicBackRadius,
			 theHistogramData -> solarBackRadius,
			 theHistogramData -> sub35SolarBackRadius);
	}
      }
      if(theHistogramData -> radius) {
#pragma omp critical
	{
	  Double_t temp[1] = {theHistogramData -> radius};
	  // ROOT is not thread safe
	  radiusTuple ->Fill(temp);
	}
      }
      delete theHistogramData;
    }// end pair for

#pragma omp critical
    {
      fractionalNuSeen += fracSeenThisLoop;
      fractionalNuThere += fracThereThisLoop;
    }
    
  }


  
  // Write out at the end
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
  
  return 0;
}
