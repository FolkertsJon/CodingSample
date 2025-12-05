// Source code for the conversion from root files out of
// the fast monte carlo to a root file containing the
// luminosity of each simulated mission 
// JF - 2/2024

#include "flightToLuminosity.hh"

namespace User
{
  
  
  // Global variable defaults
  double nOrbits = 1;
  double timeOfFlight = 0;
  double closest = 0;
  double furthest = 0;
  double kgGa = 100;
  bool ellipticalOrbit = false;
  std::string targetName = "INVALID";
  double nuPerMeVFusion = 2/26.732*kUnits::MeV; // temporarily pp1 chain again//nuPerMeVFusionFiller();
  std::vector<std::string> failedFiles = {};
  double maxRadius = 35; // controls how far out in the flight we look
  

  int singleFile(std::string filename, std::string location, TNtuple *toBeFilled){


    // get number of orbits, closest, furthest, and time step

    // Check if filename contains "EllipticalOrbit"
    ellipticalOrbit = (filename.find("EllipticalOrbit") != std::string::npos);
  
    // Replace the underscores and equal signs with spaces for easy parsing
    std::string formatted = filename;
    std::replace(formatted.begin(), formatted.end(), '_', ' ');
    std::replace(formatted.begin(), formatted.end(), '=', ' ');

    // Use a stringstream to extract the values
    std::stringstream ss(formatted);
  
    std::string tmp; // Temporary string to hold non-double values
    if(ellipticalOrbit){
      ss >> nOrbits >> tmp >> kgGa >> tmp >> tmp >> tmp >> tmp >> closest >> tmp >> furthest;
      closest = closest*kUnits::solarRadii;// add units
      furthest = furthest*kUnits::solarRadii;// add units
    }
    else{
      ss >> nOrbits >> tmp >> kgGa >> tmp >> tmp >> tmp >> tmp >> closest >> tmp >> furthest;// not right for file orbits yet
    }
    kgGa = kgGa;// add units

    // Output values for testing
    std::cout << "Num Orbits: " << nOrbits << "\n";
    std::cout << "kg gallium: " << kgGa << "\n";
    std::cout << "Closest: " << closest/kUnits::solarRadii << "\n";
    std::cout << "Furthest: " << furthest/kUnits::solarRadii << "\n";
    std::cout << "Elliptical Orbit: " << (ellipticalOrbit ? "Yes" : "No") << "\n";

  
    parseFile(filename, location, toBeFilled);
  
    return 0;
  }

  int directCalc(std::string filename, std::string location, TNtuple *toBeFilled){


    // get number of orbits, closest, furthest, and time step

    // Check if filename contains "EllipticalOrbit"
    ellipticalOrbit = (filename.find("EllipticalOrbit") != std::string::npos);
  
    // Replace the underscores and equal signs with spaces for easy parsing
    std::string formatted = filename;
    std::replace(formatted.begin(), formatted.end(), '_', ' ');
    std::replace(formatted.begin(), formatted.end(), '=', ' ');

    // Use a stringstream to extract the values
    std::stringstream ss(formatted);
  
    std::string tmp; // Temporary string to hold non-double values
    if(ellipticalOrbit){
      ss >> nOrbits >> tmp >> kgGa >> tmp >> tmp >> tmp >> tmp >> closest >> tmp >> furthest;
      closest = closest*kUnits::solarRadii;// add units
      furthest = furthest*kUnits::solarRadii;// add units
    }
    else{
      ss >> nOrbits >> tmp >> kgGa >> tmp >> tmp >> tmp >> tmp >> closest >> tmp >> furthest;// not right for file orbits yet
    }
    kgGa = kgGa;// add units

    // Output values for testing
    std::cout << "Num Orbits: " << nOrbits << "\n";
    std::cout << "kg gallium: " << kgGa << "\n";
    std::cout << "Closest: " << closest/kUnits::solarRadii << "\n";
    std::cout << "Furthest: " << furthest/kUnits::solarRadii << "\n";
    std::cout << "Elliptical Orbit: " << (ellipticalOrbit ? "Yes" : "No") << "\n";

  
    parseFileDirectCalc(filename, location, toBeFilled);
  
    return 0;
  }

  
  

  int manyFiles(std::string filePath, std::string fileNameBase, unsigned int numberOfFiles, TNtuple *toBeFilled){
    std::cout << "This mode should only be used when there is a single flight of a particular length.\n\n";
    
    
    std::cout << "What is the neutrino target: ";
    do{
      std::cin >> targetName;
      
      std::cout << "Checking if " + targetName + " is in my lists.\n";
      neutrinoTarget target = targetNameToData(targetName);
      while(target.name == "INVALID"){
	std::cout << "I don't know the name " << targetName << ". Please try another name: ";
	std::cin >> targetName;
	target = targetNameToData(targetName);
	//continue;
      }
      targetName = target.name;
       std::cout << "name = " << target.name << "\n"
		 << "PDGcode = " << target.PDGcode << "\n"
		 << "amuPerAtom = " << target.amuPerAtom << "\n"
		 << "atomicNumber = " << target.atomicNumber << "\n"
		 << "massNumber = " << target.massNumber << "\n"
		 << "nPerKilogram = " << target.nPerKilogram << "\n"
		 << "massPerAtom = " << target.massPerAtom << "\n"
		 << "targetMassFraction = " << target.targetMassFraction << "\n\n\n";
	
    }
    while(targetName == "INVALID");
  
    std::cout << "Is the orbit elliptical? (y/n): ";
    std::string temp;// for filling data
    while(temp != "y" && temp != "n"){
      std::cin >> temp;
      if(temp == "y"){
	ellipticalOrbit = true;
      }
      else if(temp == "n"){
	ellipticalOrbit = false;
      }
      else{
	std::cout << "Answer must be \"y\" or \"n\": ";
      }
    }

    if(ellipticalOrbit){
      std::cout << "What is the closest approach in Solar Radii?\n";
      std::cin >> closest;
      closest = closest*kUnits::solarRadii;// add units
    
      std::cout << "What is the furthest approach in Solar Radii? (215 =  earth, 154.8 = venus, 81.63 = mercury)\n";
      std::cin >> furthest;
      furthest = furthest*kUnits::solarRadii;// add units

      std::cout << "How many kg of target on the spacecraft?\n";
      std::cin >> kgGa;
      kgGa = 100;
    
      double a = (closest + furthest)/2;
      timeOfFlight = ( 2*M_PI*sqrt(a*a*a/( kUnits::G * kUnits::mSun ) ) );///kUnits::sec;

    }
    else{

      std::cout << "How many kg of target on the spacecraft?\n";
      std::cin >> kgGa;
      kgGa = 100;
    
    }

    std::cout << "Enter a time for the flight in months. \n(If you want the default for an ellipse or the ACO orbit, enter 0)\n";
    std:: cin >> timeOfFlight;
    timeOfFlight *= kUnits::months;
    std::cout << "Time of flight is: " << timeOfFlight/kUnits::yr << "years.\n";
  
    std::cout << "I'm parsing " << numberOfFiles << " files\n";
    double mylog = log(numberOfFiles)/log(10); // log base 10
    int modulo;
    
    double exponent = round(mylog - 2);
    if (exponent < 0) {
      // Handle the case where exponent is negative
      modulo = 1;  // This ensures no division by zero
    } else {
      modulo = (int)pow(10, exponent);
    }
    
    //#pragma omp parallel for
    std::cout << "For loop begins\n\n\n";
    for(unsigned int i = 0; i <= numberOfFiles; i++){
      //#pragma omp critical
    
      std::string filename = fileNameBase + std::to_string(i) + ".root";
      //std::cout << "The file to open is: " << filename << " at location " << filePath << "\n";
      if(i % modulo == 0){
	//std::cout << "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBH!\n";
	std::cout << "\033[F\033[FOn file " << i << " of " << numberOfFiles <<"\n";
	std::cout << "Filename = "+filename+"\n";  
      }
      
      //std::cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAH!\n";
      parseFile(filename, filePath, toBeFilled);
    }
    std::cout << "For loop ends\n";

    if(failedFiles.size() > 0){
      std::string toCout = "Several files failed to open. They were :\n";
      for(std::string file : failedFiles){
	toCout += file + "\n";
      }
      toCout = "Several files failed to open. They were ^^^^\n";
      std::cout << toCout;
    }

    return 1;
  }
  
  
  
  void parseFile(std::string filename, std::string location, TNtuple *toBeFilled){
    // std::cout << "parseFile begins.\n";
    //double amuGaPerAtom = 69.723*kUnits::amu; // number of amu in a gallium


    neutrinoTarget theTarget;
    theTarget = targetNameToData(targetName);
  
    // calculate time of flight if not already set
    if(ellipticalOrbit && timeOfFlight == 0){
      double a = (closest + furthest)/2;
      timeOfFlight = ( 2*M_PI*sqrt(a*a*a/( kUnits::G * kUnits::mSun ) ) );///kUnits::sec;
    }
    else if(timeOfFlight == 0){
      // no longer only flight; fix at some point (we rarely run in a way it's a problem)
      //timeOfFlight = 2326.33*kUnits::day;// ACO Orbit
      timeOfFlight = 3652.99809298990*kUnits::day; // Kyle Orbit
    }
  
  
    // open file
    TFile* myFile = TFile::Open((location+filename).c_str());
    if (!myFile || myFile->IsZombie()) {
      failedFiles.push_back(filename);
      //std::cout << "Could not open file " + filename + " at location " + location + "continuing.\n"
      if (myFile) {
	myFile->Close(); // Close the file only if it is not null
	//std::cout << "File exists but became a zombie.\n";
	return;
      }
      //std::cout << "File does not exist.\n";
      return;
    }
    // open tuple
  
    std::cout << "Opening Tuple\n";
    TNtuple* neutrinoTuple = new TNtuple();
    myFile->GetObject("myTuple",neutrinoTuple);
    TNtuple* radiusTuple = new TNtuple();
    myFile->GetObject("radiusTuple",radiusTuple);
  

  
    // calculate the time below maxRadius solar radii
    double timeScalar = radiusToTimeScalar(radiusTuple, maxRadius);
    if(timeScalar == -1){
      // error message is handled in the radiusToTimeScalar function.
      myFile -> Close();
      delete radiusTuple;
      delete neutrinoTuple;
      return;
    }
    std::cout << "% of time below " + std::to_string(maxRadius) + " solar radii is : " << timeScalar*100 << "\n";

    // no longer hardcoded to 35 but didn't change variable name
    double totalTimeBelow35 = timeScalar * timeOfFlight;

    double binWidth = 1e-6*kUnits::AU/kUnits::solarRadii;
    //double binWidth = 1e-1;
    double xMax = maxRadius + binWidth/2;
    double xMin = xMax;
    int nBins = 1;
    while(xMin > 0){
      nBins++;
      xMin -= binWidth;
    }
  
    // find total number of neutrinos that were there via summing over 4 pi R_event^2/total cross section
    double neutrinoTotal = 0;// need to account for fraction of neutrinos below threshold
    double uncertainty = 0;

    //double nGa = kgGa/(kUnits::amu*amuGaPerAtom)*0.339; // multiplied by percent Ga-71
    //double nPerAMU = 0.339/amuGaPerAtom;
    //double nGa = kgGa*nPerAMU; // multiplied by percent Ga-71

    double nGa = kgGa*theTarget.nPerKilogram*theTarget.targetMassFraction;

    //std::cout << "There are " << kgGa << " * " << theTarget.nPerKilogram*theTarget.targetMassFraction << " = " << nGa << " gallium atoms on this ride.\n\n";
    //std::cout << "kgGa = " << kgGa << "\n";


    // Cross sections are weighted based on the solar neutrino flux spectrum
    // we only need the weighted average
    double crossSectionPerAtom = User::bahcallEffCrossSection();

    
    double totalCrossSection = crossSectionPerAtom * nGa ;// / totalNuFlux * cutNuFlux;

    //std::cout << "The  cross section per atom in centimeters squared is : " << crossSectionPerAtom/kUnits::cm/kUnits::cm << "\n";
    //std::cout << "The number of atoms is : " << nGa << "\n";
    //std::cout << "The total cross section in meters squared is : " << totalCrossSection/kUnits::m/kUnits::m << "\n";
  
    // find total power by multiplying bt total # nus, dividing by time below maxRadius solar radii, and dividing
    // the number of neutrinos/MeV Fusion



    //std::cout << std::setprecision(15) << nuPerMeVFusion << " nu/MeV\n";
    //std::cout << std::setprecision(15) << 2/nuPerMeVFusion << " MeV/2nu\n";


    // std::cout << "The total amount of energy corresponds to " << neutrinoTotal << " / "
    //	      << nuPerMeVFusion << " = " << neutrinoTotal/nuPerMeVFusion << " MeV\n"
    //	      << nuPerMeVFusion << " = " << neutrinoTotal/nuPerMeVFusion/kUnits::J << " Joules\n";
    
    double nuPerSec = neutrinoTotal/(totalTimeBelow35/kUnits::sec);// removed hardcoded gallium frac since the orbits take care of it now.
    double nuPerSecSigma = uncertainty/(totalTimeBelow35/kUnits::sec);

    // These powers are wrong somehow, but the number of neutrinos per second is reasonable
    // double checking
  
  
    double totalPower = nuPerSec/nuPerMeVFusion;
    //std::cout << "Total Power generation is " << nuPerSec << "/" << nuPerMeVFusion << " = " << totalPower/kUnits::J*kUnits::s << " J/s\n";
    double powerUncertainty = nuPerSecSigma/nuPerMeVFusion;

    double powerInSolar = totalPower/kUnits::solarLuminosity;
    //std::cout << powerInSolar << "\%.\n";
    double percentUncertainty = powerUncertainty/totalPower;
    //std::cout << "Whhhhhhy? " << nuPerSec << "\n";
    toBeFilled -> Fill(powerInSolar, nuPerSec, percentUncertainty);

    delete neutrinoTuple;
    delete radiusTuple;
    myFile -> Close();
    // close canvas
  }



  


  void parseFileDirectCalc(std::string filename, std::string location, TNtuple *toBeFilled){

    neutrinoTarget theTarget;
    theTarget = targetNameToData(targetName);
  
      double maxRadius = 35;
    /*
    // calculate time of flight if not already set
    if(ellipticalOrbit && timeOfFlight == 0){
    double a = (closest + furthest)/2;
    timeOfFlight = ( 2*M_PI*sqrt(a*a*a/( kUnits::G * kUnits::mSun ) ) );///kUnits::sec;
    }
    else if(timeOfFlight == 0){
    // no longer only flight; fix at some point (we rarely run in a way it's a problem)
    timeOfFlight = 2326.33*kUnits::day;// only non-elliptical orbit I run currently
    }
    */
  
    // open file
    TFile* myFile = TFile::Open((location+filename).c_str());
    if (!myFile || myFile->IsZombie()) {
      myFile->Close();
      return;
    }
    // open tuple
  
    //std::cout << "Opening Tuple\n";
    TNtuple* neutrinoTuple = new TNtuple();
    myFile->GetObject("myTuple",neutrinoTuple);
  
    // open canvas
    TCanvas *c0 = new TCanvas("c0","NeutrinoHistogram",1920,1080);
    c0->cd();


  

    double binWidth = 1e-6*kUnits::AU/kUnits::solarRadii;
    //double binWidth = 1e-1;
    double xMax = maxRadius + binWidth/2;
    double xMin = xMax;
    int nBins = 1;
    while(xMin > 0){
      nBins++;
      xMin -= binWidth;
    }
  
    // find total number of neutrinos that were there via summing over 4 pi R_event^2/total cross section
    double neutrinoTotal = 0;// need to account for fraction of neutrinos below threshold
    double uncertainty = 0;


    //double nGa = kgGa/(kUnits::amu*amuGaPerAtom)*0.339; // multiplied by percent Ga-71
    //double nPerAMU = 0.339/amuGaPerAtom;
    //double nGa = kgGa*nPerAMU; // multiplied by percent Ga-71

    double nGa = kgGa*theTarget.nPerKilogram*theTarget.targetMassFraction;

    //std::cout << "There are " << kgGa << " * " << theTarget.nPerKilogram << " = " << nGa << " gallium atoms on this ride.\n\n";

    // may be depricated?
    double energyThreshold = 0;//User::galliumThreshold();

    double totalNuFlux = (User::ppFlux(0, 20*kUnits::MeV)
			  +User::pepFlux(0, 20*kUnits::MeV)
			  +User::Be7Flux(0, 20*kUnits::MeV)
			  +User::N13Flux(0, 20*kUnits::MeV)
			  +User::O15Flux(0, 20*kUnits::MeV)
			  +User::F17Flux(0, 20*kUnits::MeV)
			  +User::B8Flux(0, 20*kUnits::MeV)
			  +User::hepFlux(0, 20*kUnits::MeV)
			  );

    double cutNuFlux = (User::ppFlux(energyThreshold, 20*kUnits::MeV)
			+User::pepFlux(energyThreshold, 20*kUnits::MeV)
			+User::Be7Flux(energyThreshold, 20*kUnits::MeV)
			+User::N13Flux(energyThreshold, 20*kUnits::MeV)
			+User::O15Flux(energyThreshold, 20*kUnits::MeV)
			+User::F17Flux(energyThreshold, 20*kUnits::MeV)
			+User::B8Flux(energyThreshold, 20*kUnits::MeV)
			+User::hepFlux(energyThreshold, 20*kUnits::MeV)
			);


    // weighted average of cross section by the fusion type's flux over total flux
    // don't know what the value should actually look like. When I get the 

    // Cross sections are weighted based on the solar neutrino flux spectrum
    // we only need the weighted average
    double crossSectionPerAtom = (User::ppCrossSec(0,20*kUnits::MeV)*User::ppFlux(0, 20*kUnits::MeV)
				  + User::pepCrossSec(0,20*kUnits::MeV)*User::pepFlux(0, 20*kUnits::MeV)
				  + User::Be7CrossSec(0,20*kUnits::MeV)*User::Be7Flux(0, 20*kUnits::MeV)
				  + User::N13CrossSec(0,20*kUnits::MeV)*User::N13Flux(0, 20*kUnits::MeV)
				  + User::O15CrossSec(0,20*kUnits::MeV)*User::O15Flux(0, 20*kUnits::MeV)
				  + User::F17CrossSec(0,20*kUnits::MeV)*User::F17Flux(0, 20*kUnits::MeV)
				  + User::B8CrossSec(0,20*kUnits::MeV)*User::B8Flux(0, 20*kUnits::MeV)
				  + User::hepCrossSec(0,20*kUnits::MeV)*User::hepFlux(0, 20*kUnits::MeV))/totalNuFlux;

    
    double totalCrossSection = crossSectionPerAtom * nGa ;// / totalNuFlux * cutNuFlux;

    //std::cout << "The  cross section per atom in centimeters squared is : " << crossSectionPerAtom/kUnits::cm/kUnits::cm << "\n";
    //std::cout << "The number of atoms is : " << nGa << "\n";
    //std::cout << "The total cross section in meters squared is : " << totalCrossSection/kUnits::m/kUnits::m << "\n";

    // Uncertainty had been caculated wrong using some non-weighted stuff
    // I looked up some stats stuff, and with that I found the variance is:
    // SUM(w_i*w_i*sigma_i*sigma_i)/Sum^2(w_i) for weights w_i and uncertainties sigma_i
    // this simplifies to SS_meas/(S_meas)^2
    
    
    double totalTimeBelow35 = 0;

    // Define variables to hold the data for each entry in the TNtuple
    float timeStep, neutrinoRadius, neutrinoSignal;

    TBranch* neutrinoBranch = neutrinoTuple->GetBranch("neutrinoRadius");
    // Set the branch addresses
    neutrinoTuple->SetBranchAddress("timeStep", &timeStep);
    neutrinoTuple->SetBranchAddress("neutrinoRadius", &neutrinoRadius);
    neutrinoTuple->SetBranchAddress("neutrinoSignal", &neutrinoSignal);

    
    // Iterate over all entries in the TNtuple
    Long64_t nEntries = neutrinoTuple->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
      // Get the current entry
      neutrinoTuple->GetEntry(i);
      
      if(neutrinoRadius > maxRadius || neutrinoRadius < 1) continue;
      
      // add timeStep to totalTimeBelow35
      totalTimeBelow35 += timeStep;

      double neutrinosThisStep = neutrinoSignal * timeStep;
      
      // was dividing by the mass of a single gallium atom
      double weight = neutrinosThisStep * 4*M_PI*(neutrinoRadius*kUnits::solarRadii*neutrinoRadius*kUnits::solarRadii)/totalCrossSection;
      neutrinoTotal += weight;// sum
      // nu total uncertaintyStat squared
      uncertainty += weight*weight; // sum of squares

	
    }
    neutrinoTotal = neutrinoTotal/nOrbits; // nOrbits is 1 unless set for single files
    uncertainty = uncertainty/(neutrinoTotal * neutrinoTotal);
    //uncertainty = uncertainty/nEntries - (neutrinoTotal/nEntries)*(neutrinoTotal/nEntries);// var =  SS/N - (Sum/N)^2
    uncertainty = sqrt(uncertainty); // stDev =  sqrt(uncertainty)
    //std::cout << "The total number of neutrinos that have flown past me is: " << neutrinoTotal << " +/- "
    //	      << uncertainty << "\n\n";
  
    // find total power by multiplying bt total # nus, dividing by time below maxRadius solar radii, and dividing
    // the number of neutrinos/MeV Fusion

    std::cout << std::setprecision(15) << nuPerMeVFusion << " nu/MeV\n";
    std::cout << std::setprecision(15) << 2/nuPerMeVFusion << " MeV/2nu\n";
  
    double nuPerSec = neutrinoTotal/(totalTimeBelow35/kUnits::sec);// removed hardcoded gallium frac since the orbits take care of it now.
    double nuPerSecSigma = uncertainty/(totalTimeBelow35/kUnits::sec);

    // These powers are wrong somehow, but the number of neutrinos per second is reasonable
    // double checking
  
  
    double totalPower = nuPerSec/nuPerMeVFusion;
    double powerUncertainty = nuPerSecSigma/nuPerMeVFusion;

    double powerInSolar = totalPower/kUnits::solarLuminosity;
    //std::cout << powerInSolar << "\%.\n";
    double percentUncertainty = powerUncertainty/totalPower;

    toBeFilled -> Fill(powerInSolar, nuPerSec, percentUncertainty);

    delete neutrinoTuple;
    myFile -> Close();
    // close canvas
    delete c0;
  }









  

  
   neutrinoTarget targetNameToData(std::string theTarget){
     neutrinoTarget targetToReturn;


     if( theTarget == "Ga" ||
	 theTarget == "Nautral Gallium" ||
	 theTarget == "Gallium"){
       
       targetToReturn.name = "Nautral Gallium";        // Target's name
       targetToReturn.PDGcode = -1;         // Particle Data Group code
       targetToReturn.amuPerAtom = 69.723; // Atomic mass unit per atom
       targetToReturn.atomicNumber = 31;    // Atomic number
       targetToReturn.massNumber = 0;      // Mass number
       targetToReturn.nPerKilogram = kUnits::kg / (69.723*kUnits::amu); // number of targets per kilogram
       targetToReturn.nPerMass = kUnits::kg / (69.723*kUnits::amu) / kUnits::kg ; // number of targets per mass
       targetToReturn.massPerAtom = 69.723/kUnits::amu;  // mass per atom
       targetToReturn.targetMassFraction = 0.39892;  // mass fraction of desired target

       
     } else if (theTarget == "Ga-71" ||
		theTarget == "Ga71" ||
		theTarget == "71-Ga" ||
		theTarget == "Gallium 71" ||
		theTarget == "71Ga") {
       targetToReturn.name = "Gallium 71";        // Target's name
       targetToReturn.PDGcode = 1000310710;         // Particle Data Group code
       targetToReturn.amuPerAtom = 70.924703; // Atomic mass unit per atom
       targetToReturn.atomicNumber = 31;    // Atomic number
       targetToReturn.massNumber = 71;      // Mass number
       targetToReturn.nPerKilogram = kUnits::kg / (70.924703*kUnits::amu); // number of targets per kilogram
       targetToReturn.nPerMass = kUnits::kg / (70.924703*kUnits::amu) / kUnits::kg; // number of targets per mass
       targetToReturn.targetMassFraction = 1;  // mass fraction of desired target
       
     } else{
       targetToReturn.name = "INVALID";// Target's name
       std::cout << "That target name is invalid. I only know about Gallium and Gallium-71 currently.\n";
     }
     return targetToReturn;
   }

  double radiusToTimeScalar(TNtuple* radiusTuple, double maxRadius){
    // open canvas
    TCanvas *c0 = new TCanvas("c0","NeutrinoHistogram",1920,1080);
    c0->cd();
  
    //std::cout << "Drawing neutrinoTuple\n";
    radiusTuple->Draw("radius","radius > 1");
    auto htemp0 = (TH1F*)gPad->GetPrimitive("htemp");
    if (!htemp0) {

      // close canvas
      delete c0;
      std::cout << "NTuple was empty." << std::endl;
      return -1;
    }
    int nEntriesBig = htemp0->GetEntries();

    radiusTuple->Draw("radius",("radius > 1 && radius < " + std::to_string(maxRadius)).c_str() );
    auto htemp1 = (TH1F*)gPad->GetPrimitive("htemp");// floats are already too preciese and smaller
    int nEntriesSmall = htemp1->GetEntries();
    
    delete c0;
    
    return double(nEntriesSmall)/double(nEntriesBig);
  }




  std::pair<double,double> countNeutrinosAndUncertainty(TNtuple* neutrinoTuple, double totalCrossSection){

    Int_t nEntries = neutrinoTuple->GetEntries();

    TBranch* neutrinoBranch = neutrinoTuple->GetBranch("neutrinoRadius");
    Float_t neutrinoRadius;
    neutrinoTuple->SetBranchAddress("neutrinoRadius", &neutrinoRadius);
    double neutrinoTotal = 0;
    double uncertainty = 0;
    
    // Uncertainty had been caculated wrong using some non-weighted stuff
    // I looked up some stats stuff, and with that I found the variance is:
    // SUM(w_i*w_i*sigma_i*sigma_i)/Sum^2(w_i) for weights w_i and uncertainties sigma_i
    // this simplifies to SS_meas/(S_meas)^2

    for (Int_t i = 0; i < nEntries; ++i) {
      neutrinoTuple->GetEntry(i);
      if(neutrinoRadius > maxRadius || neutrinoRadius < 1) continue;
      if(i%1 == 0){
	//std::cout << "Radius[" << i << "] = " << neutrinoRadius << "\n";
      }
      // was dividing by the mass of a single gallium atom
      double weight =4*M_PI*(neutrinoRadius*kUnits::solarRadii)*(neutrinoRadius*kUnits::solarRadii)/totalCrossSection;
      //std::cout << "with weight = " << weight << "\n";
      neutrinoTotal += weight;// sum
      // nu total uncertaintyStat squared
      uncertainty += weight*weight; // sum of squares
    }
    neutrinoTotal = neutrinoTotal/nOrbits; // nOrbits is 1 unless set for single files
    uncertainty = uncertainty/(neutrinoTotal * neutrinoTotal);
    //uncertainty = uncertainty/nEntries - (neutrinoTotal/nEntries)*(neutrinoTotal/nEntries);// var =  SS/N - (Sum/N)^2
    uncertainty = sqrt(uncertainty); // stDev =  sqrt(uncertainty)
    //std::cout << "The total number of neutrinos that have flown past me is: " << neutrinoTotal << " +/- "
    //	      << uncertainty << "\n\n";
    return std::pair<double, double>(neutrinoTotal, uncertainty);
  }

  

  double nuPerMeVFusionFiller(){
    // This works reasonably well when I am running pp only for the cross section, but il get smaller
    // because of the higher energy fusion processes producing fewer neutrinos/MeV. I need to find a weighted
    // average of
    // I need to re-weight this based on the fractions of what neutrinos I can see. My high threshold
    // of about 400 keV cuts out a large fraction of the low-energy pp stuff, but not much of the other stuff
    double m_proton = 938.27208816*kUnits::MeV; // proton alone
    double m_electron = 510.99895*kUnits::keV;
    double m_duteron = 1.875612928*kUnits::GeV; // bare duteron
    double m_He3 = 2.80839153*kUnits::GeV; // bare helium
    double m_He4 = 3.72737933*kUnits::GeV; // bare helium
    double m_Be7 = 6.534184*kUnits::GeV; // bare berylliumn

    // Trying the by process thing
    double Epp = -(m_duteron + m_electron- 2*m_proton) -
      ( (m_He3 - m_duteron - m_proton)
	+ 0.8492*(m_He4 + 2*m_proton - 2*m_He3)
	+ 1e-7*( m_He4 + m_electron - m_proton - m_He3)
	+ 0.1508*( (m_Be7 - m_He3 - m_He4)
		   + 0.999*( 2*m_He4 - m_Be7 - m_electron - m_proton)
		   + 0.001*( 2*m_He4 - m_electron - m_Be7 - m_proton)
		   )
	);
    double Nupp = 1 +
      ( 1e-7*(1)
	+ 0.1508*(1)
	);
  
    double Epep = -(m_duteron - m_electron - 2*m_proton) -
      ( (m_He3 - m_duteron - m_proton)
	+ 0.8492*(m_He4 + 2*m_proton - 2*m_He3)
	+ 1e-7*( m_He4 + m_electron - m_proton - m_He3)
	+ 0.1508*( (m_Be7 - m_He3 - m_He4)
		   + 0.999*( 2*m_He4 - m_Be7 - m_electron - m_proton)
		   + 0.001*( 2*m_He4 - m_electron - m_Be7 - m_proton)
		   )
	);
    double Nupep = 1 +
      ( 1e-7*(1)
	+ 0.1508*(1)
	);

    double Ehep = -( m_He4 + m_electron - m_proton - m_He3) // start with the hep guaranteed energy
      -(m_He3 - m_duteron - m_proton) // include the energy from the only possible step
      -0.9975*(m_duteron + m_electron- 2*m_proton)
      -0.0025*(m_duteron - m_electron- 2*m_proton)
      ;
    double Nuhep = 1 +
      ( 0.0025*1
	+ 0.9975*(1)
	);// actually just two, but follow the logic for all of the processes
  
    double EBe7 = -( 2*m_He4 - m_Be7 - m_electron - m_proton)// only possible energy
      -(m_Be7 - m_He3 - m_He4) // more guaranteed
      -(m_He3 - m_duteron - m_proton) // more guaranteed
      -0.9975*(m_duteron + m_electron- 2*m_proton)//decision point
      -0.0025*(m_duteron - m_electron- 2*m_proton)//decision point
      ;
    double NuBe7 = 1 +
      ( 0.0025*1
	+ 0.9975*(1)
	);// actually just two, but follow the logic for all of the processes
  
    double EB8 = -( 2*m_He4 - m_electron - m_Be7 - m_proton)// only possible energy
      -(m_Be7 - m_He3 - m_He4) // more guaranteed
      -(m_He3 - m_duteron - m_proton) // more guaranteed
      -0.9975*(m_duteron + m_electron- 2*m_proton)//decision point
      -0.0025*(m_duteron - m_electron- 2*m_proton)//decision point
      ;
    double NuB8 = 1 +
      ( 0.0025*1
	+ 0.9975*(1)
	);// actually just two, but follow the logic for all of the processes
  
    double ECNO = (4*m_proton - m_He4 - 2*m_electron); // 4p -> He + 2 positron + 2 nu
    
    double NuCNO = 2;// 4p -> He + 2 positron + 2 nu

    // don't need to divide off the flux since we're dividing E by Nu. Would otherwise; did anyway
    double totalNuFlux = (User::ppFlux(0, 20*kUnits::MeV)
			  +User::pepFlux(0, 20*kUnits::MeV)
			  +User::Be7Flux(0, 20*kUnits::MeV)
			  +User::N13Flux(0, 20*kUnits::MeV)
			  +User::O15Flux(0, 20*kUnits::MeV)
			  +User::F17Flux(0, 20*kUnits::MeV)
			  +User::B8Flux(0, 20*kUnits::MeV)
			  +User::hepFlux(0, 20*kUnits::MeV)
			  );

    
    double EAll = (Epp*User::ppFlux(0, 20*kUnits::MeV)
		   +Epep*User::pepFlux(0, 20*kUnits::MeV)
		   +EBe7*User::Be7Flux(0, 20*kUnits::MeV)
		   +ECNO*User::N13Flux(0, 20*kUnits::MeV)
		   +ECNO*User::O15Flux(0, 20*kUnits::MeV)
		   +ECNO*User::F17Flux(0, 20*kUnits::MeV)
		   +EB8*User::B8Flux(0, 20*kUnits::MeV)
		   +Ehep*User::hepFlux(0, 20*kUnits::MeV)
		   )/totalNuFlux;

    double NuAll = (Nupp*User::ppFlux(0, 20*kUnits::MeV)
		    +Nupep*User::pepFlux(0, 20*kUnits::MeV)
		    +NuBe7*User::Be7Flux(0, 20*kUnits::MeV)
		    +NuCNO*User::N13Flux(0, 20*kUnits::MeV)
		    +NuCNO*User::O15Flux(0, 20*kUnits::MeV)
		    +NuCNO*User::F17Flux(0, 20*kUnits::MeV)
		    +NuB8*User::B8Flux(0, 20*kUnits::MeV)
		    +Nuhep*User::hepFlux(0, 20*kUnits::MeV)
		    )/totalNuFlux;

  

  
  
  
    return NuAll/EAll;//nuPerMeVpp*fracpp + nuPerMeVCNO*fracCNO; // ppI only, first approximation https://arxiv.org/pdf/2105.13858.pdf
  }
  
}// End User Namespace
