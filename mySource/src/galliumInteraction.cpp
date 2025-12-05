// This subprogram returns the neutrino interaction rate at 1 AU for one atom of gallium.
// This value can then be scaled up to the number of gallium atoms present, and scaled
// down with a position function. This function is designed to return a single number 
// that can be put into a constant.
// JF - 2/2021

#include "galliumInteraction.hh"

namespace User
{
  
  bool excitedOnly = true;
  bool doOscillations = false;

  std::vector<User::energyProbPair> oscillationData = User::readData("../../mySource/data/MSW_Prob_ee.txt");

  
  std::vector<double> *Fluxes = nullptr;
  
  std::vector<double>* fluxesFiller(){
    
    std::vector<double> *theReturn = new std::vector<double>;
    double eMin = 0;//galliumThreshold();
    double eMax = 20*kUnits::MeV; // 20 MeV is above the solar spectrum limit

    std::cout << "fluxesFiller() starts\n";
    
    theReturn -> push_back(hepFlux(eMin, eMax));
    theReturn -> push_back(B8Flux(eMin, eMax));
    theReturn -> push_back(F17Flux(eMin, eMax));
    theReturn -> push_back(O15Flux(eMin, eMax));
    theReturn -> push_back(N13Flux(eMin, eMax));
    theReturn -> push_back(ppFlux(eMin, eMax));
    theReturn -> push_back(pepFlux(eMin, eMax));
    theReturn -> push_back(Be7Flux(eMin, eMax));

    for(int i = 0; i < theReturn -> size(); i++){
      std::cout  << "Flux with index " << i << " = " << (*theReturn)[i] << '\n';
    }
    
    std::cout << "fluxesFiller() ends\n";
    
    
    return theReturn;// replace later
  }
  
  std::vector<double>* fluxesFiller(double eMin){
    
    std::vector<double> *theReturn = new std::vector<double>;
    double eMax = 20*kUnits::MeV;

    std::cout << "fluxesFiller(double eMin) starts\n";
    
    theReturn -> push_back(hepFlux(eMin, eMax));
    theReturn -> push_back(B8Flux(eMin, eMax));
    theReturn -> push_back(F17Flux(eMin, eMax));
    theReturn -> push_back(O15Flux(eMin, eMax));
    theReturn -> push_back(N13Flux(eMin, eMax));
    theReturn -> push_back(ppFlux(eMin, eMax));
    theReturn -> push_back(pepFlux(eMin, eMax));
    theReturn -> push_back(Be7Flux(eMin, eMax));

    for(int i = 0; i < theReturn -> size(); i++){
      std::cout  << "Flux with index " << i << " = " << (*theReturn)[i] << '\n';
    }
    
    std::cout << "fluxesFiller(double eMin) ends\n";
    
    
    return theReturn;// replace later
  }
  
  std::vector<double>* fluxReturner(){
    if(User::Fluxes == nullptr) {
      User::Fluxes = fluxesFiller();
    }
    return Fluxes;
  };
  
  std::vector<double>* fluxReturner(double eMin){
    if(User::Fluxes == nullptr) {
      User::Fluxes = fluxesFiller(eMin);
    }
    return Fluxes;
  };

  double *totalFluxPointer = nullptr;
  
  double totalFlux(){
    if( totalFluxPointer != nullptr){
      return *totalFluxPointer;
    }
    else{
      // Define a vector to hold the fluxes
      std::vector<double> fluxes = *fluxReturner();

      // Calculate the total flux
      double totalFlux = 0;
      for (double flux : fluxes) {
	totalFlux += flux;
      }
      totalFluxPointer = new double(totalFlux);
    }
    // error code
    return *totalFluxPointer;
  }
  
  
  
  // NEUTRINOS BY PROCESS AND ENERGY

  // DO NOT USE! POWER LAW IS BAD ESTIMATOR
  double powerIntegral(double xMin, double yMin, double xMax, double yMax, double intMin, double intMax) {
    // handle intMax < xMin and vice versa
    if( (intMax < xMin) || (intMin > xMax) ) return 0;
    
    // Take the natural log of the input variables
    double ln_xMin = log(xMin);
    double ln_xMax = log(xMax);
    double ln_yMin = log(yMin);
    double ln_yMax = log(yMax);
    
    // Calculate slope and intercept for the line
    double b = (ln_yMax - ln_yMin) / (ln_xMax - ln_xMin); // slope
    double ln_A = ln_yMin - b * ln_xMin; // intercept

    // restrict to the range of interest
    intMin = std::max(intMin,xMin);
    intMax = std::min(intMax,xMax);

    double integral;
    if(b+1 > 0){
      integral = exp(ln_A)*(pow(intMax,b+1)-pow(intMin,b+1))/(b+1);
    }
    else if (b+1 < 0){
      double power = -(b+1);
      integral = exp(ln_A)*(1/pow(intMax,power)-1/pow(intMin,power))/(b+1);
    }
    else{
      integral = exp(ln_A)*log(intMax/intMin);
    }
    //std::cout << "b = " << b << " eMin = " << intMin << " integral = " <<integral << "\n";
    if(isnan(integral)) integral = 0;
    return integral;
    //return LinearFitResult{slope, intercept};
  }



  
  TGraph *globalOscillationGraph = nullptr;
  
  TGraph* oscillationProbabilityGraph() {
    // Check if the graph already exists in the current ROOT session
    if (globalOscillationGraph != nullptr) {
      //std::cout << "Found the graph";
      return globalOscillationGraph;
    }

    if(oscillationData.empty()){
      oscillationData = User::readData("../../mySource/data/MSW_Prob_ee.txt");
    }
    

    // Prepare arrays to hold the energy and probability values
    const size_t nPoints = oscillationData.size();
    double *energies = new double[nPoints];
    double *probabilities = new double[nPoints];

    for (size_t i = 0; i < nPoints; ++i) {
      energies[i] = oscillationData[i].first;
      probabilities[i] = oscillationData[i].second;
      //std::cout << "E = " << energies[i] << "\nP_ee = " << probabilities[i] <<"\n";
    }//globalOscillationGraph

    // Create a new graph with the energy-probability data
    globalOscillationGraph = new TGraph(nPoints, energies, probabilities);
    globalOscillationGraph->SetName("oscillationProbabilityGraph");
    globalOscillationGraph->SetTitle("Energy vs. Probability;Energy (MeV);P_{ee}");



    // Clean up the dynamically allocated memory for arrays
    delete[] energies;
    delete[] probabilities;

    return globalOscillationGraph;
  }



  double oscillationProbability(double energy) {
    if( !doOscillations) return 1.0;
    // Retrieve the graph
    TGraph *graph = oscillationProbabilityGraph();
    
    // Use the TGraph::Eval() function to interpolate the probability for the given energy
    if (graph != nullptr) {
      //std::cout << "Oscillation prob is: " << graph->Eval(energy) << "\n";
      return graph->Eval(energy);
    } else {
      // In case the graph could not be retrieved or created, return an error value or handle appropriately
      std::cerr << "Error: Graph not available." << std::endl;
      return -1; 
    }
  }
  
  

  double hepFlux (double eMin, double eMax){
    // Found source from https://arxiv.org/abs/astro-ph/0412440 BS05(AGS,OPAL)

    double totalFlux = 8.23*1e3/kUnits::cm/kUnits::cm/kUnits::s;
    
    std::string initialPath = "../../mySource/data/hepspectrum.dat";
    std::string deeperPath = "../../../mySource/data/hepspectrum.dat";
 
    std::ifstream inFile(initialPath);
    if (!inFile) {
      inFile.open(deeperPath);
        
      if (!inFile) {
	std::cerr << "Unable to open file\n";
	exit(1); // call system to stop
      }
    }
    
    std::vector<double> E = {};
    std::vector<double> fluxes = {};
    
    std::string line;
    while (std::getline(inFile, line)) {
      std::istringstream iss(line);
      double e, flux;
      if (!(iss >> e >> flux)) { 
	break; // error in reading
      }
      if(flux == 0) flux = 1e-100;
      E.push_back(e);
      fluxes.push_back(flux);
    }
    
    inFile.close();

    double partInt = 0;
    double fullInt = 0;
    
    for(int i = 0; i < (E.size()-1); i++){
      fullInt += User::trapezoidIntegral(E[i], fluxes[i], E[i+1], fluxes[i+1], 0*kUnits::MeV, 20*kUnits::MeV);
      if(fluxes[i] == 1e-100 && fluxes[i+1]) continue; // don't include zero fluxes
      if(doOscillations){
	//std::cout << i <<"\n";
	partInt += User::trapezoidIntegral(E[i], fluxes[i]*oscillationProbability(E[i]), E[i+1], fluxes[i+1]*oscillationProbability(E[i+1]), eMin, eMax);
      }
      else{
      partInt += User::trapezoidIntegral(E[i], fluxes[i], E[i+1], fluxes[i+1], eMin, eMax);
      }
    }
    //std::cout << totalFlux*partInt/fullInt;
    return totalFlux*partInt/fullInt;
  }

  
  double B8Flux (double eMin, double eMax){
    // Found source from https://arxiv.org/abs/astro-ph/0412440 BS05(AGS,OPAL)

    double totalFlux = 4.59*1e6/kUnits::cm/kUnits::cm/kUnits::s;
    
    std::string initialPath = "../../mySource/data/b8.dat";
    std::string deeperPath = "../../../mySource/data/b8.dat";

    
    std::ifstream inFile(initialPath);
    if (!inFile) {
      inFile.open(deeperPath);
        
      if (!inFile) {
	std::cerr << "Unable to open file\n";
	exit(1); // call system to stop
      }
    }

    
    std::vector<double> E = {};
    std::vector<double> fluxes = {};
    
    std::string line;
    while (std::getline(inFile, line)) {
      std::istringstream iss(line);
      double e, flux;
      if (!(iss >> e >> flux)) { 
	break; // error in reading
      }
      E.push_back(e);
      fluxes.push_back(flux);
    }
    
    inFile.close();

    double partInt = 0;
    double fullInt = 0;
    
    for(int i = 0; i < E.size()-1; i++){
      //std::cout << partInt << ", " << E[i] <<'\n';
      fullInt += User::trapezoidIntegral(E[i], fluxes[i], E[i+1], fluxes[i+1], 0*kUnits::MeV, 20*kUnits::MeV);
      if(fluxes[i] == 1e-100 && fluxes[i+1]) continue; // don't include zero fluxes
      if(doOscillations){
	partInt += User::trapezoidIntegral(E[i], fluxes[i]*oscillationProbability(E[i]), E[i+1], fluxes[i+1]*oscillationProbability(E[i+1]), eMin, eMax);
      }
      else{
      partInt += User::trapezoidIntegral(E[i], fluxes[i], E[i+1], fluxes[i+1], eMin, eMax);
      }
    }
    //std::cout << partInt << '\n';
    
    return totalFlux*partInt/fullInt;
  }

  
  double F17Flux (double eMin, double eMax){
    // Found source from https://arxiv.org/abs/astro-ph/0412440 BS05(AGS,OPAL)

    double totalFlux = 3.31*1e6/kUnits::cm/kUnits::cm/kUnits::s;
    
    std::string initialPath = "../../mySource/data/f17.dat";
    std::string deeperPath = "../../../mySource/data/f17.dat";

    
    std::ifstream inFile(initialPath);
    if (!inFile) {
      inFile.open(deeperPath);
        
      if (!inFile) {
	std::cerr << "Unable to open file\n";
	exit(1); // call system to stop
      }
    }


    std::vector<double> E = {};
    std::vector<double> fluxes = {};
    
    std::string line;
    while (std::getline(inFile, line)) {
      std::istringstream iss(line);
      double e, flux;
      if (!(iss >> e >> flux)) { 
	break; // error in reading
      }
      E.push_back(e);
      fluxes.push_back(flux);
    }
    
    inFile.close();

    double partInt = 0;
    double fullInt = 0;
    
    for(int i = 0; i < E.size()-1; i++){
      //partInt += User::trapezoidIntegral(E[i], fluxes[i], E[i+1], fluxes[i+1], eMin, eMax);
      fullInt += User::trapezoidIntegral(E[i], fluxes[i], E[i+1], fluxes[i+1], 0*kUnits::MeV, 20*kUnits::MeV);
      if(fluxes[i] == 1e-100 && fluxes[i+1]) continue; // don't include zero fluxes
      if(doOscillations){
	partInt += User::trapezoidIntegral(E[i], fluxes[i]*oscillationProbability(E[i]), E[i+1], fluxes[i+1]*oscillationProbability(E[i+1]), eMin, eMax);
      }
      else{
      partInt += User::trapezoidIntegral(E[i], fluxes[i], E[i+1], fluxes[i+1], eMin, eMax);
      }
    }
    
    return totalFlux*partInt/fullInt;
  }

  double O15Flux (double eMin, double eMax){
    // Found source from https://arxiv.org/abs/astro-ph/0412440 BS05(AGS,OPAL)

    double totalFlux = 1.47*1e8/kUnits::cm/kUnits::cm/kUnits::s;
    
    std::string initialPath = "../../mySource/data/o15.dat";
    std::string deeperPath = "../../../mySource/data/o15.dat";

    
    std::ifstream inFile(initialPath);
    if (!inFile) {
      inFile.open(deeperPath);
        
      if (!inFile) {
	std::cerr << "Unable to open file \n";
	exit(1); // call system to stop
      }
    }

    std::vector<double> E = {};
    std::vector<double> fluxes = {};
    
    std::string line;
    while (std::getline(inFile, line)) {
      std::istringstream iss(line);
      double e, flux;
      if (!(iss >> e >> flux)) { 
	break; // error in reading
      }
      E.push_back(e);
      fluxes.push_back(flux);
    }
    
    inFile.close();

    double partInt = 0;
    double fullInt = 0;
    
    for(int i = 0; i < E.size()-1; i++){
      //partInt += User::trapezoidIntegral(E[i], fluxes[i], E[i+1], fluxes[i+1], eMin, eMax);
      fullInt += User::trapezoidIntegral(E[i], fluxes[i], E[i+1], fluxes[i+1], 0*kUnits::MeV, 20*kUnits::MeV);
      if(fluxes[i] == 1e-100 && fluxes[i+1]) continue; // don't include zero fluxes
      if(doOscillations){
	partInt += User::trapezoidIntegral(E[i], fluxes[i]*oscillationProbability(E[i]), E[i+1], fluxes[i+1]*oscillationProbability(E[i+1]), eMin, eMax);
      }
      else{
      partInt += User::trapezoidIntegral(E[i], fluxes[i], E[i+1], fluxes[i+1], eMin, eMax);
      }
    }
    
    return totalFlux*partInt/fullInt;
  }

  double N13Flux (double eMin, double eMax){
    // Found source from https://arxiv.org/abs/astro-ph/0412440 BS05(AGS,OPAL)

    double totalFlux = 2.03*1e8/kUnits::cm/kUnits::cm/kUnits::s;

    std::string initialPath = "../../mySource/data/n13.dat";
    std::string deeperPath = "../../../mySource/data/n13.dat";

    
    std::ifstream inFile(initialPath);
    if (!inFile) {
      inFile.open(deeperPath);
        
      if (!inFile) {
	std::cerr << "Unable to open file \n";
	exit(1); // call system to stop
      }
    }

    std::vector<double> E = {};
    std::vector<double> fluxes = {};
    
    std::string line;
    while (std::getline(inFile, line)) {
      std::istringstream iss(line);
      double e, flux;
      if (!(iss >> e >> flux)) { 
	break; // error in reading
      }
      E.push_back(e);
      fluxes.push_back(flux);
    }
    
    inFile.close();

    double partInt = 0;
    double fullInt = 0;
    
    for(int i = 0; i < E.size()-1; i++){
      //partInt += User::trapezoidIntegral(E[i], fluxes[i], E[i+1], fluxes[i+1], eMin, eMax);
      fullInt += User::trapezoidIntegral(E[i], fluxes[i], E[i+1], fluxes[i+1], 0*kUnits::MeV, 20*kUnits::MeV);
      if(fluxes[i] == 1e-100 && fluxes[i+1]) continue; // don't include zero fluxes
      if(doOscillations){
	partInt += User::trapezoidIntegral(E[i], fluxes[i]*oscillationProbability(E[i]), E[i+1], fluxes[i+1]*oscillationProbability(E[i+1]), eMin, eMax);
      }
      else{
      partInt += User::trapezoidIntegral(E[i], fluxes[i], E[i+1], fluxes[i+1], eMin, eMax);
      }
    }
    
    return totalFlux*partInt/fullInt;
  }

  double ppFlux (double eMin, double eMax){
    // Found source from https://arxiv.org/abs/astro-ph/0412440 BS05(AGS,OPAL)

    double totalFlux =  6.05*1e10/kUnits::cm/kUnits::cm/kUnits::s;
    
    std::string initialPath = "../../mySource/data/ppenergytab.txt";
    std::string deeperPath = "../../../mySource/data/ppenergytab.txt";

    
    std::ifstream inFile(initialPath);
    if (!inFile) {
      inFile.open(deeperPath);
        
      if (!inFile) {
	std::cerr << "Unable to open file\n";
	exit(1); // call system to stop
      }
    }

    std::vector<std::pair<double, double>> data;

    std::string line;
    while (std::getline(inFile, line)) {
      std::istringstream iss(line);
      double e, flux;
      while (iss >> e >> flux) {
        if (iss.fail()) {
	  iss.clear();
	  std::string tmp;
	  iss >> tmp; // read and discard non-numeric string
        }
        else {
	  data.push_back(std::make_pair(e, flux));
        }
      }
    }
    
    inFile.close();

    // Sort the data vector based on the first element of the pair
    std::sort(data.begin(), data.end());

    // Now, split the sorted data back into separate E and fluxes vectors
    std::vector<double> E, fluxes;
    for(const auto& pair : data) {
      E.push_back(pair.first);
      fluxes.push_back(pair.second);
    }

    double partInt = 0;
    double fullInt = 0;
    
    for(int i = 0; i < E.size()-1; i++){
      //partInt += User::trapezoidIntegral(E[i], fluxes[i], E[i+1], fluxes[i+1], eMin, eMax);
      fullInt += User::trapezoidIntegral(E[i], fluxes[i], E[i+1], fluxes[i+1], 0*kUnits::MeV, 20*kUnits::MeV);
      if(fluxes[i] == 1e-100 && fluxes[i+1]) continue; // don't include zero fluxes
      if(doOscillations){
	partInt += User::trapezoidIntegral(E[i], fluxes[i]*oscillationProbability(E[i]), E[i+1], fluxes[i+1]*oscillationProbability(E[i+1]), eMin, eMax);
      }
      else{
	partInt += User::trapezoidIntegral(E[i], fluxes[i], E[i+1], fluxes[i+1], eMin, eMax);
      }
    }
    
    return totalFlux*partInt/fullInt;
  }  

  double pepFlux(double eMin, double eMax){
    // Found source from https://arxiv.org/abs/astro-ph/0412440 BS05(AGS,OPAL)
    // Doesn't currently account for fraction that's below threshold
    if(eMin < 1.442234*kUnits::MeV && eMax > 1.442234*kUnits::MeV){
      return 1.45*1e8/kUnits::cm/kUnits::cm/kUnits::s
	*oscillationProbability(1.442234*kUnits::MeV);
    }
    else{
      return 0;
    }
  }

  double Be7Flux(double eMin, double eMax){
    // Found source from https://arxiv.org/abs/astro-ph/0412440 BS05(AGS,OPAL)
    // Doesn't currently account for fraction that's below threshold

    // Check if maximum energy is less than or equal to minimum energy
    if(eMax <= eMin) {
        throw std::invalid_argument("emin must be less than emax");
    }

    if(eMin > 0.861*kUnits::MeV || eMax < 0.383*kUnits::MeV){ // neither line
      return 0;
    }
    double flux = 4.38*1e9/kUnits::cm/kUnits::cm/kUnits::s;
    if(eMin > 0.383*kUnits::MeV && eMax >= 0.861){ // Lower energy line
      return (flux - 5e8/kUnits::cm/kUnits::cm/kUnits::s)
	*oscillationProbability(0.861*kUnits::MeV);
    }
    if(eMin <= 0.383*kUnits::MeV && eMax < 0.861){// Higher energy line
      return 5e8/kUnits::cm/kUnits::cm/kUnits::s
	*oscillationProbability(0.383*kUnits::MeV);
    }
    double oscillationCorrection =(
				   (5e8/kUnits::cm/kUnits::cm/kUnits::s)
				   *oscillationProbability(0.383*kUnits::MeV)
				   +(flux-5e8/kUnits::cm/kUnits::cm/kUnits::s)
				   *oscillationProbability(0.861*kUnits::MeV)
				   )/flux;
    return flux * oscillationCorrection;// both lines
  }

  // found source to make integration not matter
  
 
  
  // Cross sections for processes

  double ppCrossSec(double eMin, double eMax){
    //std::cout << "excitedOnly: " << excitedOnly << "\n";
    if(excitedOnly){
      return 0.00493044e-46*kUnits::cm*kUnits::cm; // from bahcall code ignoring ground
    }
    return 11.6104e-46*kUnits::cm*kUnits::cm; // from bahcall code
  }

  double pepCrossSec(double eMin, double eMax){
    if(excitedOnly){
      return 37.275e-46*kUnits::cm*kUnits::cm; // from bahcall code ignoring ground
    }
    return 203.902e-46*kUnits::cm*kUnits::cm; // from bahcall code
  }

  // I should split into 7Be1 and 7Be2 for the two lines like in Bahcall code
  double Be7CrossSec(double eMin, double eMax){
    if(eMin <= 390*kUnits::keV){
      if(eMin <= 233.2*kUnits::keV && !excitedOnly){ // below both means full cross sec (check exact number, only approximate rn
	return 71.776e-46*kUnits::cm*kUnits::cm;// from bahcall code averaged by flux weighting
      }
      return 4.079e-46*kUnits::cm*kUnits::cm; // from bahcall code ignoring ground  averaged by flux weighting
    }
    return 0;
  }

  double N13CrossSec(double eMin, double eMax){
    if(excitedOnly){
      return 3.62366e-46*kUnits::cm*kUnits::cm; // from bahcall code ignoring ground
    }
    return 60.4193e-46*kUnits::cm*kUnits::cm; // from bahcall code averaged by flux weighting
    //return 60.4e-46/kUnits::cm/kUnits::cm; // from first approx
  }

  double O15CrossSec(double eMin, double eMax){
    if(excitedOnly){
      return 15.8765e-46*kUnits::cm*kUnits::cm; // from bahcall code ignoring ground
    }
    return 113.696e-46*kUnits::cm*kUnits::cm; // from bahcall code averaged by flux weighting
  }

  double F17CrossSec(double eMin, double eMax){
    if(excitedOnly){
      return 16.0581e-46*kUnits::cm*kUnits::cm; // from bahcall code ignoring ground
    }
    return 114.365e-46*kUnits::cm*kUnits::cm; // from bahcall code averaged by flux weighting
    //return 113.9e-46/kUnits::cm/kUnits::cm; // from first approx
  }

  double B8CrossSec(double eMin, double eMax){
    if(excitedOnly){
      return 21222.9e-46*kUnits::cm*kUnits::cm; // from bahcall code ignoring ground
    }
    return 24028.7e-46*kUnits::cm*kUnits::cm; // from bahcall code averaged by flux weighting
    //return 2.4e-42/kUnits::cm/kUnits::cm; // from first approx
  }

  double hepCrossSec(double eMin, double eMax){
    if(excitedOnly){
      return 66110.7e-46*kUnits::cm*kUnits::cm; // from bahcall code ignoring ground
    }
    return 71431.5e-46*kUnits::cm*kUnits::cm; // from bahcall code averaged by flux weighting
  }


  double bahcallEffCrossSection(){
    double totalNuFlux = (User::ppFlux(0, 20*kUnits::MeV)
			  +User::pepFlux(0, 20*kUnits::MeV)
			  +User::Be7Flux(0, 20*kUnits::MeV)
			  +User::N13Flux(0, 20*kUnits::MeV)
			  +User::O15Flux(0, 20*kUnits::MeV)
			  +User::F17Flux(0, 20*kUnits::MeV)
			  +User::B8Flux(0, 20*kUnits::MeV)
			  +User::hepFlux(0, 20*kUnits::MeV)
			  );
    return (User::ppCrossSec(0,20*kUnits::MeV)*User::ppFlux(0, 20*kUnits::MeV)
      + User::pepCrossSec(0,20*kUnits::MeV)*User::pepFlux(0, 20*kUnits::MeV)
      + User::Be7CrossSec(0,20*kUnits::MeV)*User::Be7Flux(0, 20*kUnits::MeV)
      + User::N13CrossSec(0,20*kUnits::MeV)*User::N13Flux(0, 20*kUnits::MeV)
      + User::O15CrossSec(0,20*kUnits::MeV)*User::O15Flux(0, 20*kUnits::MeV)
      + User::F17CrossSec(0,20*kUnits::MeV)*User::F17Flux(0, 20*kUnits::MeV)
      + User::B8CrossSec(0,20*kUnits::MeV)*User::B8Flux(0, 20*kUnits::MeV)
      + User::hepCrossSec(0,20*kUnits::MeV)*User::hepFlux(0, 20*kUnits::MeV))/totalNuFlux;
  }


  // Generate random Energy for processes
  // Define a structure to hold energy and flux
  struct EnergyFlux {
    double energy;
    double flux;
  };

  // Define a comparison function for sorting
  bool compareEnergy(const EnergyFlux& a, const EnergyFlux& b) {
    return a.energy < b.energy;
  }


  double ppEnergy(TH1D* hist){// need to get actual number
    // Create a random number generator
    TRandom3 *randGen = new TRandom3(0);

    double sampledEnergy = 0;
    while(sampledEnergy < galliumThreshold()){
      // Randomly sample an energy from the histogram
      sampledEnergy = hist->GetRandom(randGen);
    }

    // Clean up
    delete randGen;

    return sampledEnergy;
  }

  
  TGraph* ppEnergyGraphMaker(){
    

    // Open the data file
    std::string initialPath = "../../mySource/data/ppenergytab.txt";
    std::string deeperPath = "../../../mySource/data/ppenergytab.txt";

    
    std::ifstream inFile(initialPath);
    if (!inFile) {
      inFile.open(deeperPath);
        
      if (!inFile) {
	std::cerr << "Unable to open file \n";
	exit(1); // call system to stop
      }
    }

    // Create a vector for energy-flux pairs
    std::vector<EnergyFlux> energyFluxPairs;

    // Read the data from the file
    std::string line;
    while (std::getline(inFile, line)) {
      if (line.empty() || line[0] == ' ' || line[0] == '\t' || line[0] == 'q') {
	continue;
      }
      std::stringstream ss(line);
      double e, f;
      while (ss >> e >> f) {
	EnergyFlux ef = {e * kUnits::MeV, f};
	energyFluxPairs.push_back(ef);
      }
    }

    // Sort the energyFluxPairs based on energy
    std::sort(energyFluxPairs.begin(), energyFluxPairs.end(), compareEnergy);

    // Extract sorted energy and flux vectors
    std::vector<double> sortedEnergy;
    std::vector<double> sortedFlux;

    // exponential extension down to 100 eV for all graphs
    double x0 = energyFluxPairs[0].energy;
    double x1 = energyFluxPairs[1].energy;
    double y0 = energyFluxPairs[0].flux;
    double y1 = energyFluxPairs[1].flux;
    double b = log(y0/y1)/(x0-x1);
    double logA = log(y0)-b*x0;

    // pp goes the smallest so skip
    //sortedEnergy.push_back(100*kUnits::eV);
    //sortedFlux.push_back( exp(logA) * exp(b*100*kUnits::eV) );
	
    for (const auto& ef : energyFluxPairs) {
      sortedEnergy.push_back(ef.energy);
      sortedFlux.push_back(ef.flux);
    }
    inFile.close();
      
    double normalization = 0;
      
    for(int i = 0; i < sortedEnergy.size(); i++){
      if(i == 0) continue; // skip first for integration
	
	
      normalization += trapezoidIntegral(sortedEnergy[i-1], sortedFlux[i-1], sortedEnergy[i], sortedFlux[i], 0*kUnits::MeV, 20*kUnits::MeV);
    }
      

    for(int i = 0; i < sortedFlux.size(); i++){
      sortedFlux[i] *= ppFlux(0*kUnits::MeV, 20*kUnits::MeV)/normalization*kUnits::cm*kUnits::cm*kUnits::s;
    }

    // Create arrays from the vectors
    const int numPoints = sortedEnergy.size();
    double* xArray = &sortedEnergy[0];
    double* yArray = &sortedFlux[0];
      
    // Create the graph
    TGraph *tg = new TGraph(numPoints, xArray, yArray);
    tg->SetTitle("PP Flux Spectrum");
    tg->GetXaxis()->SetTitle("Energy (MeV)");
    tg->GetYaxis()->SetTitle("Flux (cm^{-2} s^{-1} MeV^{-1})");
    tg->SetLineColor(kYellow);
    tg->SetMarkerColor(kYellow);
      
      

     
    auto c1 = new TCanvas("ppCanvas","canvas1",1200,900);
    c1->cd();
    c1->SetLogy(1);
    c1->SetLogx(1);

    tg->Draw();
    c1->Print("ppSpectrum.pdf");
    
    delete c1;
    
    
    return tg;
  }
    

  TH1D* ppEnergyHistogramMaker(){
    TH1D* hist = (TH1D*)gDirectory->Get("ppHist");
    
    // If the histogram already exists, return it
    if (hist != NULL) return hist;
    {
      //#pragma omp critical
      //std::cout << "in critical pp";
      //hist = (TH1D*)gDirectory->Get("ppHist");
      //if (hist != NULL) return hist;
    

      // Open the data file
    std::string initialPath = "../../mySource/data/ppenergytab.txt";
    std::string deeperPath = "../../../mySource/data/ppenergytab.txt";

    
    std::ifstream inFile(initialPath);
    if (!inFile) {
      std::cerr << "Unable to open file at: " << initialPath << ". Trying one level deeper...\n";
      inFile.open(deeperPath);
        
      if (!inFile) {
	std::cerr << "Unable to open file at: " << deeperPath << "\n";
	exit(1); // call system to stop
      }
    }

      // Create a vector for energy-flux pairs
      std::vector<EnergyFlux> energyFluxPairs;

      // Read the data from the file
      std::string line;
      while (std::getline(inFile, line)) {
	if (line.empty() || line[0] == ' ' || line[0] == '\t' || line[0] == 'q') {
	  continue;
	}
	std::stringstream ss(line);
	double e, f;
	while (ss >> e >> f) {
	  EnergyFlux ef = {e * kUnits::MeV, f};
	  energyFluxPairs.push_back(ef);
	  //std::cout << "Energy = " << e << " Flux = " << f <<'\n';
	}
      }

      // Sort the energyFluxPairs based on energy
      std::sort(energyFluxPairs.begin(), energyFluxPairs.end(), compareEnergy);

      // Extract sorted energy and flux vectors
      std::vector<double> sortedEnergy;
      std::vector<double> sortedFlux;

      // exponential extension down to 100 eV for all graphs
      double x0 = energyFluxPairs[0].energy;
      double x1 = energyFluxPairs[1].energy;
      double y0 = energyFluxPairs[0].flux;
      double y1 = energyFluxPairs[1].flux;
      double b = log(y0/y1)/(x0-x1);
      double logA = log(y0)-b*x0;
      
      sortedEnergy.push_back(100*kUnits::eV);
      sortedFlux.push_back( exp(logA) * exp(b*100*kUnits::eV) );
      for (const auto& ef : energyFluxPairs) {
	sortedEnergy.push_back(ef.energy);
	sortedFlux.push_back(ef.flux);
      }
      inFile.close();

      // Create a histogram
      hist = new TH1D("ppHist", "Energy Spectrum", 24, sortedEnergy.back(), sortedEnergy.back());//

      // Fill the histogram with fluxes
      for (unsigned int i = 0; i < sortedEnergy.size(); i++) {
	hist->Fill(sortedEnergy[i],sortedFlux[i]);
      }
      auto c1 = new TCanvas("ppCanvas","canvas1",1200,900);
      c1->cd();
      c1->SetLogy(1);
      c1->SetLogx(1);

      hist->Draw();
      c1->Print("ppSpectrum.pdf");
    
      delete c1;
    
    }
    return hist;
  }

  
  double pepEnergy(){// need to get actual number
    return 1.442*kUnits::MeV; 
  }

  
  double Be7Energy(){// should be close enough
    double fracLowEnergy = 4e8/(4.38e9);
    
    // Generate a random number in the range [0,1)
    TRandom3 randGenerator(0); // generate from system clock
    double energy = 0;
    while(energy < galliumThreshold()){
      double random = randGenerator.Uniform(1);

      if( random < fracLowEnergy){
	energy = 383*kUnits::keV;
      }
      else{
	energy = 861*kUnits::keV;
      }
    }
    return energy;
  }
  
  TGraph* N13EnergyGraphMaker(){
    
      //if (hist != NULL) return hist;
    
      // Open the data file
      std::ifstream inFile("../../mySource/data/n13.dat");
      if (!inFile) {
	std::cerr << "Unable to open file ppenergytab.txt\n";
	exit(1);   // call system to stop
      }

      // Create a vector for energy-flux pairs
      std::vector<EnergyFlux> energyFluxPairs;

      // Read the data from the file
      std::string line;
      while (std::getline(inFile, line)) {
	std::stringstream ss(line);
	double e, f;
	while (ss >> e >> f) {
	  EnergyFlux ef = {e * kUnits::MeV, f};
	  energyFluxPairs.push_back(ef);
	}
      }

      // Sort the energyFluxPairs based on energy
      std::sort(energyFluxPairs.begin(), energyFluxPairs.end(), compareEnergy);

      // Extract sorted energy and flux vectors
      std::vector<double> sortedEnergy;
      std::vector<double> sortedFlux;

      // exponential extension down to 100 eV for all graphs
      double x0 = energyFluxPairs[0].energy;
      double x1 = energyFluxPairs[1].energy;
      double y0 = energyFluxPairs[0].flux;
      double y1 = energyFluxPairs[1].flux;
      double b = log(y0/y1)/(x0-x1);
      double logA = log(y0)-b*x0;
      
      sortedEnergy.push_back(100*kUnits::eV);
      sortedFlux.push_back( exp(logA) * exp(b*100*kUnits::eV) );
      
      
      for (const auto& ef : energyFluxPairs) {
	sortedEnergy.push_back(ef.energy);
	sortedFlux.push_back(ef.flux);
      }
      inFile.close();
      
      double normalization = 0;
      
      for(int i = 0; i < sortedEnergy.size(); i++){
	if(i == 0) continue; // skip first for integration
	
	
	normalization += trapezoidIntegral(sortedEnergy[i-1], sortedFlux[i-1], sortedEnergy[i], sortedFlux[i], 0*kUnits::MeV, 20*kUnits::MeV);
      }
      

      for(int i = 0; i < sortedFlux.size(); i++){
	sortedFlux[i] *= N13Flux(0*kUnits::MeV, 20*kUnits::MeV)/normalization*kUnits::cm*kUnits::cm*kUnits::s;
      }

      // Create arrays from the vectors
      const int numPoints = sortedEnergy.size();
      double* xArray = &sortedEnergy[0];
      double* yArray = &sortedFlux[0];
      
      // Create the graph
      TGraph *tg = new TGraph(numPoints, xArray, yArray);
      tg->SetTitle("N13 Flux Spectrum");
      tg->GetXaxis()->SetTitle("Energy (MeV)");
      tg->GetYaxis()->SetTitle("Flux (cm^{-2} s^{-1} MeV^{-1})");
      tg->SetLineColor(kRed);
      tg->SetMarkerColor(kRed);
      
      

     
      auto c1 = new TCanvas("N13Canvas","canvas1",1200,900);
      c1->cd();
      c1->SetLogy(1);
      c1->SetLogx(1);

      tg->Draw("A L P *");
      c1->Print("N13Spectrum.pdf");
    
      delete c1;
    
    
    return tg;
  }
  
  TH1D* N13EnergyHistogramMaker(){
    TH1D* hist = (TH1D*)gDirectory->Get("N13Hist");
    
    // If the histogram already exists, return it
    if (hist != NULL) return hist;
    
      //if (hist != NULL) return hist;
    
      // Open the data file
      std::ifstream inFile("../../mySource/data/n13.dat");
      if (!inFile) {
	std::cerr << "Unable to open file ppenergytab.txt\n";
	exit(1);   // call system to stop
      }

      // Create a vector for energy-flux pairs
      std::vector<EnergyFlux> energyFluxPairs;

      // Read the data from the file
      std::string line;
      while (std::getline(inFile, line)) {
	std::stringstream ss(line);
	double e, f;
	while (ss >> e >> f) {
	  EnergyFlux ef = {e * kUnits::MeV, f};
	  energyFluxPairs.push_back(ef);
	}
      }

      // Sort the energyFluxPairs based on energy
      std::sort(energyFluxPairs.begin(), energyFluxPairs.end(), compareEnergy);

      // Extract sorted energy and flux vectors
      std::vector<double> sortedEnergy;
      std::vector<double> sortedFlux;

      // exponential extension down to 100 eV for all graphs
      double x0 = energyFluxPairs[0].energy;
      double x1 = energyFluxPairs[1].energy;
      double y0 = energyFluxPairs[0].flux;
      double y1 = energyFluxPairs[1].flux;
      double b = log(y0/y1)/(x0-x1);
      double logA = log(y0)-b*x0;
      
      sortedEnergy.push_back(100*kUnits::eV);
      sortedFlux.push_back( exp(logA) * exp(b*100*kUnits::eV) );
      
      for (const auto& ef : energyFluxPairs) {
	sortedEnergy.push_back(ef.energy);
	sortedFlux.push_back(ef.flux);
      }
      inFile.close();
      {	
	//#pragma omp critical
      // Create a histogram
      hist = new TH1D("N13Hist", "Energy Spectrum", 24, sortedEnergy.front(), sortedEnergy.back());

      // Fill the histogram with fluxes
      for (unsigned int i = 0; i < sortedEnergy.size(); i++) {
	hist->Fill(sortedEnergy[i],sortedFlux[i]);
      }
      auto c1 = new TCanvas("N13Canvas","canvas1",1200,900);
      c1->cd();
      c1->SetLogy(1);
      c1->SetLogx(1);

      hist->Draw();
      c1->Print("N13Spectrum.pdf");
    
      delete c1;

    }
    return hist;
  }
  
  TGraph* O15EnergyGraphMaker(){
    
      // Open the data file
      std::ifstream inFile("../../mySource/data/o15.dat");
      if (!inFile) {
	std::cerr << "Unable to open file ppenergytab.txt\n";
	exit(1);   // call system to stop
      }

      // Create a vector for energy-flux pairs
      std::vector<EnergyFlux> energyFluxPairs;

      // Read the data from the file
      std::string line;
      while (std::getline(inFile, line)) {
	std::stringstream ss(line);
	double e, f;
	while (ss >> e >> f) {
	  EnergyFlux ef = {e * kUnits::MeV, f};
	  energyFluxPairs.push_back(ef);
	}
      }

      // Sort the energyFluxPairs based on energy
      std::sort(energyFluxPairs.begin(), energyFluxPairs.end(), compareEnergy);

      // Extract sorted energy and flux vectors
      std::vector<double> sortedEnergy;
      std::vector<double> sortedFlux;

      // exponential extension down to 100 eV for all graphs
      double x0 = energyFluxPairs[0].energy;
      double x1 = energyFluxPairs[1].energy;
      double y0 = energyFluxPairs[0].flux;
      double y1 = energyFluxPairs[1].flux;
      double b = log(y0/y1)/(x0-x1);
      double logA = log(y0)-b*x0;
      
      sortedEnergy.push_back(100*kUnits::eV);
      sortedFlux.push_back( exp(logA) * exp(b*100*kUnits::eV) );
      
      for (const auto& ef : energyFluxPairs) {
	sortedEnergy.push_back(ef.energy);
	sortedFlux.push_back(ef.flux);
      }
      inFile.close();
      
      double normalization = 0;
      
      for(int i = 0; i < sortedEnergy.size(); i++){
	if(i == 0) continue; // skip first for integration
	
	
	normalization += trapezoidIntegral(sortedEnergy[i-1], sortedFlux[i-1], sortedEnergy[i], sortedFlux[i], 0*kUnits::MeV, 20*kUnits::MeV);
      }
      

      for(int i = 0; i < sortedFlux.size(); i++){
	sortedFlux[i] *= O15Flux(0*kUnits::MeV, 20*kUnits::MeV)/normalization*kUnits::cm*kUnits::cm*kUnits::s;
      }

      // Create arrays from the vectors
      const int numPoints = sortedEnergy.size();
      double* xArray = &sortedEnergy[0];
      double* yArray = &sortedFlux[0];
      
      // Create the graph
      TGraph *tg = new TGraph(numPoints, xArray, yArray);
      tg->SetTitle("O15 Flux Spectrum");
      tg->GetXaxis()->SetTitle("Energy (MeV)");
      tg->GetYaxis()->SetTitle("Flux (cm^{-2} s^{-1} MeV^{-1})");
      tg->SetLineColor(kGreen);
      tg->SetMarkerColor(kGreen);
      
      

     
      auto c1 = new TCanvas("O15Canvas","canvas1",1200,900);
      c1->cd();
      c1->SetLogy(1);
      c1->SetLogx(1);

      tg->Draw("A L P *");
      c1->Print("O15Spectrum.pdf");
    
      delete c1;
    
    
    return tg;
  }
    
  
  TH1D* O15EnergyHistogramMaker(){
    TH1D* hist = (TH1D*)gDirectory->Get("O15Hist");
    
    // If the histogram already exists, return it
    if (hist != NULL) return hist;
      //if (hist != NULL) return hist;
    
      // Open the data file
      std::ifstream inFile("../../mySource/data/o15.dat");
      if (!inFile) {
	std::cerr << "Unable to open file ppenergytab.txt\n";
	exit(1);   // call system to stop
      }

      // Create a vector for energy-flux pairs
      std::vector<EnergyFlux> energyFluxPairs;

      // Read the data from the file
      std::string line;
      while (std::getline(inFile, line)) {
	std::stringstream ss(line);
	double e, f;
	while (ss >> e >> f) {
	  EnergyFlux ef = {e * kUnits::MeV, f};
	  energyFluxPairs.push_back(ef);
	}
      }

      // Sort the energyFluxPairs based on energy
      std::sort(energyFluxPairs.begin(), energyFluxPairs.end(), compareEnergy);

      // Extract sorted energy and flux vectors
      std::vector<double> sortedEnergy;
      std::vector<double> sortedFlux;

      // exponential extension down to 100 eV for all graphs
      double x0 = energyFluxPairs[0].energy;
      double x1 = energyFluxPairs[1].energy;
      double y0 = energyFluxPairs[0].flux;
      double y1 = energyFluxPairs[1].flux;
      double b = log(y0/y1)/(x0-x1);
      double logA = log(y0)-b*x0;
      
      sortedEnergy.push_back(100*kUnits::eV);
      sortedFlux.push_back( exp(logA) * exp(b*100*kUnits::eV) );
      
      for (const auto& ef : energyFluxPairs) {
	sortedEnergy.push_back(ef.energy);
	sortedFlux.push_back(ef.flux);
      }
      inFile.close();

    {
      //#pragma omp critical
      // Create a histogram
      hist = new TH1D("O15Hist", "Energy Spectrum", 24, sortedEnergy.back(), sortedEnergy.front());

      // Fill the histogram with fluxes
      for (unsigned int i = 0; i < sortedEnergy.size(); i++) {
	hist->Fill(sortedEnergy[i],sortedFlux[i]);
      }
      auto c1 = new TCanvas("O15Canvas","canvas1",1200,900);
      c1->cd();
      c1->SetLogy(1);
      c1->SetLogx(1);

      hist->Draw();
      c1->Print("O15Spectrum.pdf");
    
      delete c1;

    }
    return hist;
  }

  TGraph* F17EnergyGraphMaker(){
    
      // Open the data file
      std::ifstream inFile("../../mySource/data/f17.dat");
      if (!inFile) {
	std::cerr << "Unable to open file ppenergytab.txt\n";
	exit(1);   // call system to stop
      }

      // Create a vector for energy-flux pairs
      std::vector<EnergyFlux> energyFluxPairs;

      // Read the data from the file
      std::string line;
      while (std::getline(inFile, line)) {
	std::stringstream ss(line);
	double e, f;
	while (ss >> e >> f) {
	  EnergyFlux ef = {e * kUnits::MeV, f};
	  energyFluxPairs.push_back(ef);
	}
      }

      // Sort the energyFluxPairs based on energy
      std::sort(energyFluxPairs.begin(), energyFluxPairs.end(), compareEnergy);

      // Extract sorted energy and flux vectors
      std::vector<double> sortedEnergy;
      std::vector<double> sortedFlux;

      // exponential extension down to 100 eV for all graphs
      double x0 = energyFluxPairs[0].energy;
      double x1 = energyFluxPairs[1].energy;
      double y0 = energyFluxPairs[0].flux;
      double y1 = energyFluxPairs[1].flux;
      double b = log(y0/y1)/(x0-x1);
      double logA = log(y0)-b*x0;
      
      sortedEnergy.push_back(100*kUnits::eV);
      sortedFlux.push_back( exp(logA) * exp(b*100*kUnits::eV) );
      
      for (const auto& ef : energyFluxPairs) {
	sortedEnergy.push_back(ef.energy);
	sortedFlux.push_back(ef.flux);
      }
      inFile.close();
      
      double normalization = 0;
      
      for(int i = 0; i < sortedEnergy.size(); i++){
	if(i == 0) continue; // skip first for integration
	
	
	normalization += trapezoidIntegral(sortedEnergy[i-1], sortedFlux[i-1], sortedEnergy[i], sortedFlux[i], 0*kUnits::MeV, 20*kUnits::MeV);
      }
      

      for(int i = 0; i < sortedFlux.size(); i++){
	sortedFlux[i] *= F17Flux(0*kUnits::MeV, 20*kUnits::MeV)/normalization*kUnits::cm*kUnits::cm*kUnits::s;
      }

      // Create arrays from the vectors
      const int numPoints = sortedEnergy.size();
      double* xArray = &sortedEnergy[0];
      double* yArray = &sortedFlux[0];
      
      // Create the graph
      TGraph *tg = new TGraph(numPoints, xArray, yArray);
      tg->SetTitle("F17 Flux Spectrum");
      tg->GetXaxis()->SetTitle("Energy (MeV)");
      tg->GetYaxis()->SetTitle("Flux (cm^{-2} s^{-1} MeV^{-1})");
      tg->SetLineColor(kAzure-6);
      tg->SetMarkerColor(kAzure-6);
      
      

     
      auto c1 = new TCanvas("O15Canvas","canvas1",1200,900);
      c1->cd();
      c1->SetLogy(1);
      c1->SetLogx(1);

      tg->Draw("A L P *");
      c1->Print("F17Spectrum.pdf");
    
      delete c1;
    
    
    return tg;
  }
    

  TH1D* F17EnergyHistogramMaker(){
    TH1D* hist = (TH1D*)gDirectory->Get("F17Hist");
    
    // If the histogram already exists, return it
    if (hist != NULL) return hist;
      //if (hist != NULL) return hist;
    
      // Open the data file
      std::ifstream inFile("../../mySource/data/f17.dat");
      if (!inFile) {
	std::cerr << "Unable to open file ppenergytab.txt\n";
	exit(1);   // call system to stop
      }

      // Create a vector for energy-flux pairs
      std::vector<EnergyFlux> energyFluxPairs;

      // Read the data from the file
      std::string line;
      while (std::getline(inFile, line)) {
	std::stringstream ss(line);
	double e, f;
	while (ss >> e >> f) {
	  EnergyFlux ef = {e * kUnits::MeV, f};
	  energyFluxPairs.push_back(ef);
	}
      }

      // Sort the energyFluxPairs based on energy
      std::sort(energyFluxPairs.begin(), energyFluxPairs.end(), compareEnergy);

      // Extract sorted energy and flux vectors
      std::vector<double> sortedEnergy;
      std::vector<double> sortedFlux;

      // exponential extension down to 100 eV for all graphs
      double x0 = energyFluxPairs[0].energy;
      double x1 = energyFluxPairs[1].energy;
      double y0 = energyFluxPairs[0].flux;
      double y1 = energyFluxPairs[1].flux;
      double b = log(y0/y1)/(x0-x1);
      double logA = log(y0)-b*x0;
      
      sortedEnergy.push_back(100*kUnits::eV);
      sortedFlux.push_back( exp(logA) * exp(b*100*kUnits::eV) );
      
      for (const auto& ef : energyFluxPairs) {
	sortedEnergy.push_back(ef.energy);
	sortedFlux.push_back(ef.flux);
      }
      inFile.close();

    {
      //#pragma omp critical
      // Create a histogram
      hist = new TH1D("F17Hist", "Energy Spectrum", 24, sortedEnergy.front(), sortedEnergy.back());

      // Fill the histogram with fluxes
      for (unsigned int i = 0; i < sortedEnergy.size(); i++) {
	hist->Fill(sortedEnergy[i],sortedFlux[i]);
      }
      auto c1 = new TCanvas("F17Canvas","canvas1",1200,900);
      c1->cd();
      c1->SetLogy(1);
      c1->SetLogx(1);

      hist->Draw();
      c1->Print("F17Spectrum.pdf");
    
      delete c1;

    }
    return hist;
  }
  
  TGraph* B8EnergyGraphMaker(){
      // Open the data file
      std::ifstream inFile("../../mySource/data/b8.dat");
      if (!inFile) {
	std::cerr << "Unable to open file ppenergytab.txt\n";
	exit(1);   // call system to stop
      }

      // Create a vector for energy-flux pairs
      std::vector<EnergyFlux> energyFluxPairs;

      // Read the data from the file
      std::string line;
      while (std::getline(inFile, line)) {
	std::stringstream ss(line);
	double e, f;
	while (ss >> e >> f) {
	  EnergyFlux ef = {e * kUnits::MeV, f};
	  energyFluxPairs.push_back(ef);
	}
      }

      // Sort the energyFluxPairs based on energy
      std::sort(energyFluxPairs.begin(), energyFluxPairs.end(), compareEnergy);

      // Extract sorted energy and flux vectors
      std::vector<double> sortedEnergy;
      std::vector<double> sortedFlux;

      // exponential extension down to 100 eV for all graphs
      double x0 = energyFluxPairs[1].energy;
      double x1 = energyFluxPairs[2].energy;
      double y0 = energyFluxPairs[1].flux;
      double y1 = energyFluxPairs[2].flux;
      double b = log(y0/y1)/(x0-x1);
      double logA = log(y0)-b*x0;
      //std::cout << "Flux = " << exp(logA) << " * " << "exp(" << b << " * x)\n"; 
      
      //sortedEnergy.push_back(100*kUnits::eV);
      //sortedFlux.push_back( exp(logA) * exp(b*100*kUnits::eV) );
      
     
      for (const auto& ef : energyFluxPairs) {
	sortedEnergy.push_back(ef.energy);
	sortedFlux.push_back(ef.flux);
      }
      // replace poinnt 0,0 with my exponential
      sortedEnergy[0] = 100*kUnits::eV;
      sortedFlux[0]= exp(logA) * exp(b*100*kUnits::eV);
          
      inFile.close();
      
      double normalization = 0;
      
      for(int i = 0; i < sortedEnergy.size(); i++){
	if(i == 0) continue; // skip first for integration
	
	
	normalization += trapezoidIntegral(sortedEnergy[i-1], sortedFlux[i-1], sortedEnergy[i], sortedFlux[i], 0*kUnits::MeV, 20*kUnits::MeV);
      }
      

      for(int i = 0; i < sortedFlux.size(); i++){
	sortedFlux[i] *= B8Flux(0*kUnits::MeV, 20*kUnits::MeV)/normalization*kUnits::cm*kUnits::cm*kUnits::s;
      }

      // Create arrays from the vectors
      const int numPoints = sortedEnergy.size();
      double* xArray = &sortedEnergy[0];
      double* yArray = &sortedFlux[0];
      
      // Create the graph
      TGraph *tg = new TGraph(numPoints, xArray, yArray);
      tg->SetTitle("B8 Flux Spectrum");
      tg->GetXaxis()->SetTitle("Energy (MeV)");
      tg->GetYaxis()->SetTitle("Flux (cm^{-2} s^{-1} MeV^{-1})");
      tg->SetLineColor(kMagenta);
      tg->SetMarkerColor(kMagenta);
      
      

     
      auto c1 = new TCanvas("B8Canvas","canvas1",1200,900);
      c1->cd();
      c1->SetLogy(1);
      c1->SetLogx(1);

      tg->Draw("A L P *");
      c1->Print("B8Spectrum.pdf");
    
      delete c1;
    
    
    return tg;
  }
  
  TH1D* B8EnergyHistogramMaker(){
    TH1D* hist = (TH1D*)gDirectory->Get("B8Hist");
    
    // If the histogram already exists, return it
    if (hist != NULL) return hist;
      //if (hist != NULL) return hist;
    
      // Open the data file
      std::ifstream inFile("../../mySource/data/b8.dat");
      if (!inFile) {
	std::cerr << "Unable to open file ppenergytab.txt\n";
	exit(1);   // call system to stop
      }

      // Create a vector for energy-flux pairs
      std::vector<EnergyFlux> energyFluxPairs;

      // Read the data from the file
      std::string line;
      while (std::getline(inFile, line)) {
	std::stringstream ss(line);
	double e, f;
	while (ss >> e >> f) {
	  EnergyFlux ef = {e * kUnits::MeV, f};
	  energyFluxPairs.push_back(ef);
	}
      }

      // Sort the energyFluxPairs based on energy
      std::sort(energyFluxPairs.begin(), energyFluxPairs.end(), compareEnergy);

      // Extract sorted energy and flux vectors
      std::vector<double> sortedEnergy;
      std::vector<double> sortedFlux;

      // exponential extension down to 100 eV for all graphs
      double x0 = energyFluxPairs[0].energy;
      double x1 = energyFluxPairs[1].energy;
      double y0 = energyFluxPairs[0].flux;
      double y1 = energyFluxPairs[1].flux;
      double b = log(y0/y1)/(x0-x1);
      double logA = log(y0)-b*x0;
      
      sortedEnergy.push_back(100*kUnits::eV);
      sortedFlux.push_back( exp(logA) * exp(b*100*kUnits::eV) );
      
      for (const auto& ef : energyFluxPairs) {
	sortedEnergy.push_back(ef.energy);
	sortedFlux.push_back(ef.flux);
      }
      inFile.close();

    {
      //#pragma omp critical
      // Create a histogram
      hist = new TH1D("B8Hist", "Energy Spectrum", 24, sortedEnergy.front(), sortedEnergy.back());

      // Fill the histogram with fluxes
      for (unsigned int i = 0; i < sortedEnergy.size(); i++) {
	hist->Fill(sortedEnergy[i],sortedFlux[i]);
      }
      auto c1 = new TCanvas("B8Canvas","canvas1",1200,900);
      c1->cd();
      c1->SetLogy(1);
      c1->SetLogx(1);

      hist->Draw();
      c1->Print("B8Spectrum.pdf");
    
      delete c1;

    }
    return hist;
  }

  TGraph* hepEnergyGraphMaker(){
      // Open the data file
      std::ifstream inFile("../../mySource/data/hepspectrum.dat");
      if (!inFile) {
	std::cerr << "Unable to open file hepspectrum.dat\n";
	exit(1);   // call system to stop
      }

      // Create a vector for energy-flux pairs
      std::vector<EnergyFlux> energyFluxPairs;

      // Read the data from the file
      std::string line;
      while (std::getline(inFile, line)) {
	std::stringstream ss(line);
	double e, f;
	while (ss >> e >> f) {
	  EnergyFlux ef = {e * kUnits::MeV, f};
	  energyFluxPairs.push_back(ef);
	}
      }

      // Sort the energyFluxPairs based on energy
      std::sort(energyFluxPairs.begin(), energyFluxPairs.end(), compareEnergy);

      // Extract sorted energy and flux vectors
      std::vector<double> sortedEnergy;
      std::vector<double> sortedFlux;

      // exponential extension down to 100 eV for all graphs
      double x0 = energyFluxPairs[0].energy;
      double x1 = energyFluxPairs[1].energy;
      double y0 = energyFluxPairs[0].flux;
      double y1 = energyFluxPairs[1].flux;
      double b = log(y0/y1)/(x0-x1);
      double logA = log(y0)-b*x0;
      
      sortedEnergy.push_back(100*kUnits::eV);
      sortedFlux.push_back( exp(logA) * exp(b*100*kUnits::eV) );
      
      for (const auto& ef : energyFluxPairs) {
	sortedEnergy.push_back(ef.energy);
	sortedFlux.push_back(ef.flux);
      }
      inFile.close();
      
      double normalization = 0;
      
      for(int i = 0; i < sortedEnergy.size(); i++){
	if(i == 0) continue; // skip first for integration
	
	
	normalization += trapezoidIntegral(sortedEnergy[i-1], sortedFlux[i-1], sortedEnergy[i], sortedFlux[i], 0*kUnits::MeV, 20*kUnits::MeV);
      }
      

      for(int i = 0; i < sortedFlux.size(); i++){
	sortedFlux[i] *= hepFlux(0*kUnits::MeV, 20*kUnits::MeV)/normalization*kUnits::cm*kUnits::cm*kUnits::s;
      }

      // Create arrays from the vectors
      const int numPoints = sortedEnergy.size();
      double* xArray = &sortedEnergy[0];
      double* yArray = &sortedFlux[0];
      
      // Create the graph
      TGraph *tg = new TGraph(numPoints, xArray, yArray);
      tg->SetTitle("hep Flux Spectrum");
      tg->GetXaxis()->SetTitle("Energy (MeV)");
      tg->GetYaxis()->SetTitle("Flux (cm^{-2} s^{-1} MeV^{-1})");
      tg->SetLineColor(kCyan);
      tg->SetMarkerColor(kCyan);
      
      

     
      auto c1 = new TCanvas("hepCanvas","canvas1",1200,900);
      c1->cd();
      c1->SetLogy(1);
      c1->SetLogx(1);

      tg->Draw("A L P *");
      c1->Print("hepSpectrum.pdf");
    
      delete c1;
    
    
    return tg;
  }
    

  TH1D* hepEnergyHistogramMaker() {
    //std::cout << "Getting hepHist\n";
    TH1D* hist = (TH1D*)gDirectory->Get("hepHist");
    //std::cout << "got hepHist\n";
    
    // If the histogram already exists, return it
    if (hist != NULL) return hist;
    std::cout << "need new hepHist\n";
      //if (hist != NULL) return hist;
    
      // Open the data file
      std::ifstream inFile("../../mySource/data/hepspectrum.dat");
      if (!inFile) {
	std::cerr << "Unable to open file ppenergytab.txt\n";
	exit(1);   // call system to stop
      }

      // Create a vector for energy-flux pairs
      std::vector<EnergyFlux> energyFluxPairs;

      // Read the data from the file
      std::string line;
      while (std::getline(inFile, line)) {
	std::stringstream ss(line);
	double e, f;
	while (ss >> e >> f) {
	  EnergyFlux ef = {e * kUnits::MeV, f};
	  energyFluxPairs.push_back(ef);
	}
      }

      // Sort the energyFluxPairs based on energy
      std::sort(energyFluxPairs.begin(), energyFluxPairs.end(), compareEnergy);

      // Extract sorted energy and flux vectors
      std::vector<double> sortedEnergy;
      std::vector<double> sortedFlux;

      // exponential extension down to 100 eV for all graphs
      double x0 = energyFluxPairs[0].energy;
      double x1 = energyFluxPairs[1].energy;
      double y0 = energyFluxPairs[0].flux;
      double y1 = energyFluxPairs[1].flux;
      double b = log(y0/y1)/(x0-x1);
      double logA = log(y0)-b*x0;
      
      sortedEnergy.push_back(100*kUnits::eV);
      sortedFlux.push_back( exp(logA) * exp(b*100*kUnits::eV) );
      
      for (const auto& ef : energyFluxPairs) {
	sortedEnergy.push_back(ef.energy);
	sortedFlux.push_back(ef.flux);
      }
      inFile.close();

    {
      //#pragma omp critical
      // Create a histogram
      hist = new TH1D("hepHist", "Energy Spectrum", 24, sortedEnergy.front(), sortedEnergy.back());

      // Fill the histogram with fluxes
      for (unsigned int i = 0; i < sortedEnergy.size(); i++) {
	hist->Fill(sortedEnergy[i],sortedFlux[i]);
      }
      auto c1 = new TCanvas("hepCanvas","canvas1",1200,900);
      c1->cd();
      c1->SetLogy(1);
      c1->SetLogx(1);

      hist->Draw("H");
      c1->Print("hepSpectrum.pdf");
    
      delete c1;

    }
    return hist;
  }
  
  TGraph* pepEnergyGraphMaker(){
    // Parameters for the Gaussian
    double amplitude = pepFlux(0*kUnits::MeV, 20*kUnits::MeV)*kUnits::cm*kUnits::cm*kUnits::s;
    double mean = 1.442234*kUnits::MeV;
    double sigma = mean*1e-3;

    // Create the Gaussian function
    TF1* gaussian = new TF1("gaussian", "[0]*TMath::Gaus(x, [1], [2], kTRUE)", mean - 5 * sigma, mean + 5 * sigma); // kTRUE uses normalized gaus
    gaussian->SetParameters(amplitude, mean, sigma);

    // Number of points
    int nPoints = 103;

    // Create a TGraph
    TGraph* tg = new TGraph(nPoints);
    tg->SetTitle("pep Flux Spectrum");
    tg->GetXaxis()->SetTitle("Energy (MeV)");
    tg->GetYaxis()->SetTitle("Flux (cm^{-2} s^{-1} MeV^{-1})");
    tg->SetLineColor(kViolet-5);
    tg->SetMarkerColor(kViolet-5);

    // Calculate the step size
    double stepSize = (10 * sigma) / (nPoints - 1);

    // Evaluate the function at evenly spaced points
    //std::cout << "Evaluating pep function\n";
    for (int i = 0; i < nPoints; i++) {
        double x = mean - 5 * sigma + i * stepSize;
        double y = gaussian->Eval(x);
	if(i == 0 || i == (nPoints - 1) ) y = 0;// go to zero at extremes, loses almost nothing
	//std::cout << "(" << x << ", " << y << ") Press return to continue";
	//std::cin.ignore();
        tg->SetPoint(i, x, y);
    }
    
      auto c1 = new TCanvas("pepCanvas","canvas1",1200,900);
      c1->cd();
      c1->SetLogy(1);
      c1->SetLogx(1);

      tg->Draw("A L P *");
      c1->Print("pepSpectrum.pdf");
    
      delete c1;

    return tg;
  }

  TGraph* Be7EnergyGraphMaker(){
    // Parameters for the Gaussians
    double amplitude1 = Be7Flux(0.5*kUnits::MeV, 20*kUnits::MeV)*kUnits::cm*kUnits::cm*kUnits::s;
    double mean1 = 0.861*kUnits::MeV;
    double sigma1 = mean1*1e-3;//mean1*1e-3;
    
    double amplitude0 = Be7Flux(0*kUnits::MeV, 0.5*kUnits::MeV)*kUnits::cm*kUnits::cm*kUnits::s;
    double mean0 = 0.383*kUnits::MeV;
    double sigma0 = mean0*1e-3;

    // Create the Gaussian functions 
    TF1* gaussian0 = new TF1("gaussian0", "[0]*TMath::Gaus(x, [1], [2], kTRUE)", mean0 - 5 * sigma0, mean0 + 5 * sigma0); // kTRUE uses normalized gaussian
    gaussian0->SetParameters(amplitude0, mean0, sigma0);
    
    TF1* gaussian1 = new TF1("gaussian1", "[0]*TMath::Gaus(x, [1], [2], kTRUE)", mean1 - 5 * sigma1, mean1 + 5 * sigma1);
    gaussian1->SetParameters(amplitude1, mean1, sigma1);

    // Number of points per gaussian
    int nPoints = 103;

    // Create a TGraph
    TGraph* tg0 = new TGraph(nPoints);
    TGraph* tg1 = new TGraph(nPoints);

    // Calculate the step size
    double stepSize0 = (10 * sigma0) / (nPoints - 1);
    double stepSize1 = (10 * sigma1) / (nPoints - 1);

    // Evaluate the function at evenly spaced points
    for (int i = 0; i < nPoints; i++) {
      double x = mean0 - 5 * sigma0 + i * stepSize0;
      double y = gaussian0->Eval(x);
      if(i == 0 || i == (nPoints - 1) ) y = 0;// go to zero at extremes, loses almost nothing
      tg0->SetPoint(i, x, y);
      x = mean1 - 5 * sigma1 + i * stepSize1;
      y = gaussian1->Eval(x);
      if(i == 0 || i == (nPoints - 1) ) y = 0;// go to zero at extremes, loses almost nothing
      tg1->SetPoint(i, x, y);
    }
    
    TGraph* tg = new TGraph(tg0->GetN() + tg1->GetN());
    tg->SetTitle("Be7 Flux Spectrum");
    tg->GetXaxis()->SetTitle("Energy (MeV)");
    tg->GetYaxis()->SetTitle("Flux (cm^{-2} s^{-1} MeV^{-1})");
    tg->SetLineColor(kOrange);
    tg->SetMarkerColor(kOrange);

    // Add points from tg1
    for (int i = 0; i < tg1->GetN(); i++) {
      double x, y;
      tg1->GetPoint(i, x, y);
      tg->SetPoint(i + tg0->GetN(), x, y);
    }

    // Add points from tg0
    for (int i = 0; i < tg0->GetN(); i++) {
      double x, y;
      tg0->GetPoint(i, x, y);
      tg->SetPoint(i, x, y);
    }
    
    
    
    auto c1 = new TCanvas("Be7Canvas","canvas1",1200,900);
    c1->cd();
    c1->SetLogy(1);
    c1->SetLogx(1);

    tg->Draw("A L P *");
    c1->Print("Be7Spectrum.pdf");
    
    delete c1;

    return tg;
  }
  
  double N13Energy(TH1D* hist){// need to get actual number
    // Create a random number generator
    TRandom3 *randGen = new TRandom3(0);

    double sampledEnergy = 0;
    while(sampledEnergy < galliumThreshold()){
      // Randomly sample an energy from the histogram
      sampledEnergy = hist->GetRandom(randGen);
    }

    // Clean up
    delete randGen;

    return sampledEnergy;
  }

  double O15Energy(TH1D* hist){// need to get actual number
    // Create a random number generator
    TRandom3 *randGen = new TRandom3(0);

    double sampledEnergy = 0;
    while(sampledEnergy < galliumThreshold()){
      // Randomly sample an energy from the histogram
      sampledEnergy = hist->GetRandom(randGen);
    }

    // Clean up
    delete randGen;

    return sampledEnergy;
  }

  double F17Energy(TH1D* hist){// need to get actual number
    // Create a random number generator
    TRandom3 *randGen = new TRandom3(0);

    double sampledEnergy = 0;
    while(sampledEnergy < galliumThreshold()){
      // Randomly sample an energy from the histogram
      sampledEnergy = hist->GetRandom(randGen);
    }

    // Clean up
    delete randGen;

    return sampledEnergy;
  }

  double B8Energy(TH1D* hist){// need to get actual number
    // Create a random number generator
    TRandom3 *randGen = new TRandom3(0);

    double sampledEnergy = 0;
    while(sampledEnergy < galliumThreshold()){
      // Randomly sample an energy from the histogram
      sampledEnergy = hist->GetRandom(randGen);
    }

    // Clean up
    delete randGen;

    return sampledEnergy;
  }

  double hepEnergy(TH1D* hist){// need to get actual number
    // Create a random number generator
    TRandom3 *randGen = new TRandom3(0);

    double sampledEnergy = 0;
    while(sampledEnergy < galliumThreshold()){
      // Randomly sample an energy from the histogram
      sampledEnergy = hist->GetRandom(randGen);
    }

    // Clean up
    delete randGen;

    return sampledEnergy;
  }

  // Returns the per atom neutrino rate @ earth
  double neutrinoRate(){
    double totalFlux = 0;
    for(double flux : *fluxReturner()){
      totalFlux += flux;
    }
    std::cout << "The total flux of the neutrino processes is: "
	      << totalFlux*kUnits::cm*kUnits::cm*kUnits::s << " cm-2 s-1\n";
    std::cout << "The effective solar neutrino cross section for this flux is: "
	      << bahcallEffCrossSection()/kUnits::cm/kUnits::cm << " cm2\n";
    
    double rate = bahcallEffCrossSection()*totalFlux;
    std::cout << "With a total effective per-atom rate of :" << rate/kUnits::s
	      << " per second.\n";
    
    return rate;
  }

  // for backward compatibility
  // Returns the per atom neutrino rate @ earth
  double neutrinoRate(double eMin, double eMax){
    return neutrinoRateByProcess(eMin, eMax);
  }

  double neutrinoRateByProcess(double eMin, double eMax){
    // number of interactions per second per atom at 1 AU. 
    // Calculated from effective cross sections and solar
    // neutrino flux by source graph. 
    double neutrinoInteractionRate = 0;
    double previous = 0;
    // Integral rate
    std::cout << "At the Earth: \n";
    neutrinoInteractionRate += ppCrossSec(eMin, eMax) * ppFlux(eMin, eMax);
    std::cout << "pp rate is : " << ppCrossSec(eMin, eMax) 
	      << " * " << ppFlux(eMin, eMax) << " = "
	      << neutrinoInteractionRate - previous
	      << " neutrinos per second.\n";
    previous = neutrinoInteractionRate;
    neutrinoInteractionRate += pepCrossSec(eMin, eMax) * pepFlux(eMin, eMax);
    std::cout << "pep rate is : " << neutrinoInteractionRate - previous
	      << " neeutrinos per second.\n";
    previous = neutrinoInteractionRate;
    neutrinoInteractionRate += Be7CrossSec(eMin, eMax) * Be7Flux(eMin, eMax);
    std::cout << "7Be rate is : " << neutrinoInteractionRate - previous
	      << " neutrinos per second.\n";
    previous = neutrinoInteractionRate;
    neutrinoInteractionRate += N13CrossSec(eMin, eMax) * N13Flux(eMin, eMax);
    std::cout << "13N rate is : " << neutrinoInteractionRate - previous
	      << " neutrinos per second.\n";
    previous = neutrinoInteractionRate;
    neutrinoInteractionRate += O15CrossSec(eMin, eMax) * O15Flux(eMin, eMax);
    std::cout << "15O rate is : " << neutrinoInteractionRate - previous
	      << " neutrinos per second.\n";
    previous = neutrinoInteractionRate;
    neutrinoInteractionRate += F17CrossSec(eMin, eMax) * F17Flux(eMin, eMax);
    std::cout << "17F rate is : " << neutrinoInteractionRate - previous
	      << " neutrinos per second.\n";
    previous = neutrinoInteractionRate;
    neutrinoInteractionRate += B8CrossSec(eMin, eMax) * B8Flux(eMin, eMax);
    std::cout << "8B rate is : " << neutrinoInteractionRate - previous
	      << " neutrinos per second.\n";
    previous = neutrinoInteractionRate;
    neutrinoInteractionRate += hepCrossSec(eMin, eMax) * hepFlux(eMin, eMax);
    std::cout << "hep rate is : " << neutrinoInteractionRate - previous
	      << " neutrinos per second.\n";

    std::cout << "Total rate is : " << neutrinoInteractionRate
	      << " neutrinos per second.\n\n";
    
  
    return neutrinoInteractionRate;
  }


  // This function approximates the amplification of nearer positions having higher probabilities of electrons.
  // We assume that the factor is 2 at 1 rSol, and go linearly up from 1 at 35 rSol
  double linOscFactor(double Radius){
    if(Radius/kUnits::solarRadii > 35){
      return 1;// further than 35 solar radii
    }
    else if (Radius/kUnits::solarRadii > 1){
      double m = -1/34;
      double b = 2 + 1/34;
      double x = Radius/kUnits::solarRadii;
      double y = m*x+b;
      return y;//between 1 and 25 solar radii
    }
    return 2;//less than 1 solar radius for complete defintion.
  }
  

  
  
  // random energy for neutrino
  double neutrinoEnergy(){
    // Define a vector to hold the fluxes
    std::vector<double> fluxes = *fluxReturner();


    double eMin = galliumThreshold();
    double eMax = 20*kUnits::MeV;

    // Calculate the total flux
    double totalFlux = 0;
    for (double flux : fluxes) {
      totalFlux += flux;
    }

    // Generate a random number in the range [0, totalFlux)
    TRandom3 randGenerator(0); // generate from system clock
    double randomFlux = randGenerator.Uniform(totalFlux);

    // Find which flux the random number corresponds to
    double accumulatedFlux = 0;
    for (int i = 0; i < fluxes.size(); i++) {
      //std::cout << fluxes[i] << ',';
      accumulatedFlux += fluxes[i];
      if (randomFlux < accumulatedFlux) {
	// This is the selected flux
	switch(i) {
	case 0:
	  //std::cout << "neutrino was hep\n";
	  return hepEnergy(hepEnergyHistogramMaker());
	  break;
	case 1:
	  //std::cout << "neutrino was b8\n";
	  return B8Energy(B8EnergyHistogramMaker());
	  break;
	case 2:
	  //std::cout << "neutrino was F17\n";
	  return F17Energy(F17EnergyHistogramMaker());	  
	  break;
	case 3:
	  //std::cout << "neutrino was O15\n";
	  return O15Energy(O15EnergyHistogramMaker());	
	  break;
	case 4:
	  //std::cout << "neutrino was N13\n";
	  return N13Energy(N13EnergyHistogramMaker());
	  break;
	case 5:
	  //std::cout << "neutrino was pp\n";
	  return ppEnergy(ppEnergyHistogramMaker());
	  break;
	case 6:
	  //std::cout << "neutrino was pep\n";
	  return pepEnergy();
	  break;
	case 7:
	  //std::cout << "neutrino was Be7\n";
	  return Be7Energy();
	  break;
	default:
	  std::cout << "i is not between 0 and 8" << std::endl;
	}
	
      }
    }
    std::cout << '\n';

    return 0;  // Should not reach here
  } 






  void checkGalliumCrossSection(){

    bool temp = excitedOnly;
    bool tempOsc = doOscillations;
    
      
    
    for(bool oscil : {true, false}){
      doOscillations = oscil;
      std::cout << "doOscillations = " << doOscillations << "\n\n";
	
      double totalNuFlux = 0;
    
      totalNuFlux += User::ppFlux(0, 20*kUnits::MeV);
      std::cout << "ppFlux: " << User::ppFlux(0, 20*kUnits::MeV)*kUnits::cm*kUnits::cm*kUnits::s << " cm-2 s-1\n";
    
      totalNuFlux += User::pepFlux(0, 20*kUnits::MeV);
      std::cout << "pepFlux: " << User::pepFlux(0, 20*kUnits::MeV)*kUnits::cm*kUnits::cm*kUnits::s << " cm-2 s-1\n";
    
      totalNuFlux += User::Be7Flux(0, 20*kUnits::MeV);
      std::cout << "Be7Flux: " << User::Be7Flux(0, 20*kUnits::MeV)*kUnits::cm*kUnits::cm*kUnits::s << " cm-2 s-1\n";
    
      totalNuFlux += User::N13Flux(0, 20*kUnits::MeV);
      std::cout << "N13Flux: " << User::N13Flux(0, 20*kUnits::MeV)*kUnits::cm*kUnits::cm*kUnits::s << " cm-2 s-1\n";
    
      totalNuFlux += User::O15Flux(0, 20*kUnits::MeV);
      std::cout << "O15Flux: " << User::O15Flux(0, 20*kUnits::MeV)*kUnits::cm*kUnits::cm*kUnits::s << " cm-2 s-1\n";
    
      totalNuFlux += User::F17Flux(0, 20*kUnits::MeV);
      std::cout << "F17Flux: " << User::F17Flux(0, 20*kUnits::MeV)*kUnits::cm*kUnits::cm*kUnits::s << " cm-2 s-1\n";
    
      totalNuFlux += User::B8Flux(0, 20*kUnits::MeV);
      std::cout << "B8Flux: " << User::B8Flux(0, 20*kUnits::MeV)*kUnits::cm*kUnits::cm*kUnits::s << " cm-2 s-1\n";
    
      totalNuFlux += User::hepFlux(0, 20*kUnits::MeV);
      std::cout << "hepFlux: " << User::hepFlux(0, 20*kUnits::MeV)*kUnits::cm*kUnits::cm*kUnits::s << " cm-2 s-1\n";

    
      std::cout << "\nTotal neutrino Flux: " << totalNuFlux*kUnits::cm*kUnits::cm*kUnits::s << " cm-2 s-1\n\n";
    



      for(bool excited : {true, false}){
	excitedOnly = excited;
	std::cout << "excitedOnly = " << excitedOnly << "\n\n";
    
	double weightedCrossSection = 0;
    
	weightedCrossSection += User::ppCrossSec(0,20*kUnits::MeV)*User::ppFlux(0, 20*kUnits::MeV);
	std::cout << "ppCrossSec: " << User::ppCrossSec(0,20*kUnits::MeV)/kUnits::cm/kUnits::cm << " cm2\n";
    
	weightedCrossSection += User::pepCrossSec(0,20*kUnits::MeV)*User::pepFlux(0, 20*kUnits::MeV);
	std::cout << "pepCrossSec: " << User::pepCrossSec(0,20*kUnits::MeV)/kUnits::cm/kUnits::cm << " cm2\n";
    
	weightedCrossSection += User::Be7CrossSec(0,20*kUnits::MeV)*User::Be7Flux(0, 20*kUnits::MeV);
	std::cout << "Be7CrossSec: " << User::Be7CrossSec(0,20*kUnits::MeV)/kUnits::cm/kUnits::cm << " cm2\n";
    
	weightedCrossSection += User::N13CrossSec(0,20*kUnits::MeV)*User::N13Flux(0, 20*kUnits::MeV);
	std::cout << "N13CrossSec: " << User::N13CrossSec(0,20*kUnits::MeV)/kUnits::cm/kUnits::cm << " cm2\n";
    
	weightedCrossSection += User::O15CrossSec(0,20*kUnits::MeV)*User::O15Flux(0, 20*kUnits::MeV);
	std::cout << "O15CrossSec: " << User::O15CrossSec(0,20*kUnits::MeV)/kUnits::cm/kUnits::cm << " cm2\n";
    
	weightedCrossSection += User::F17CrossSec(0,20*kUnits::MeV)*User::F17Flux(0, 20*kUnits::MeV);
	std::cout << "F17CrossSec: " << User::F17CrossSec(0,20*kUnits::MeV)/kUnits::cm/kUnits::cm << " cm2\n";
    
	weightedCrossSection += User::B8CrossSec(0,20*kUnits::MeV)*User::B8Flux(0, 20*kUnits::MeV);
	std::cout << "B8CrossSec: " << User::B8CrossSec(0,20*kUnits::MeV)/kUnits::cm/kUnits::cm << " cm2\n";
    
	weightedCrossSection += User::hepCrossSec(0,20*kUnits::MeV)*User::hepFlux(0, 20*kUnits::MeV);
	std::cout << "hepCrossSec: " << User::hepCrossSec(0,20*kUnits::MeV)/kUnits::cm/kUnits::cm << " cm2\n";

	weightedCrossSection /= totalNuFlux;
	std::cout << "\nWeighted Cross section: " << weightedCrossSection/kUnits::cm/kUnits::cm << " cm2\n\n";
    
      }
    }

    excitedOnly = temp;
    doOscillations = tempOsc;
    
  } // end void checking gallium


  
}// End User Namespace
