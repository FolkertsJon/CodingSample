// This file contains all the functions needed to call from main.cc
// in order to change the parameters (time, filemode, etc)
// This contains the functions that can check for arguments, parse
// the argument(s) that are input, and then execute the appropriate
// changes.
// Jonathan Folkerts: Created May 2021
// 5/12/2021 : exists; added option parsing 
// 5/14/2021 : added Par struct for holding paramters and functions 
//             to initialize it and change values in it

#include <variableParser.hh>

namespace User{

  bool printedHelp = false;
  
  void parInitializer(User::Par *p){
    p->fileMode = false;
    //p->myFile = "../ACOorbits.txt";
    //std::cout << p->myFile << " is the name of the file.\n";
    p->timeStep= 1*kUnits::min;// time step for iterating
    p->detectorRadius = 0.3*kUnits::meter; // approximate value
    p->cosmicAcceptance = 0.00001; // 5 nines
    p->solarAcceptance = p->cosmicAcceptance;// same as cosmic
    p->kgGallium = 25*kUnits::kg; //25 kg default
    p->galliumAtomNumber = p->kgGallium * 8.637e24/kUnits::kg;// calculate gallium atom number
    p->closest = 10*kUnits::solarRadii; // Closest approach for elliptical mode
    p->furthest = 1.1*kUnits::venusSunDistance; // Furthest approach for elliptical mode
    p->nLoops = 100; // number of iterations to run
    p->solarNeutrinoRate = 0;//User::neutrinoRate(0,20)*p->galliumAtomNumber*0.399; // per atom rate * nAtoms * percentage gallium 71
    p->cosBack = User::cosmicBackground(100,100); // cosmic background rate
    p->inclination = 1; // currently unused
    p->timeLimit = 0;
    p->excitedOnly = 0;
  }

  void parChanger(User::Par *p, User::Variable v){
    
    if(v.varName == "fileMode") {
      p->fileMode = (v.varVal);
      std::cout << "fileMode is " << p->fileMode << std::endl;
      /*if(p->fileMode){
	std::string theDefault = "../ACOorbits.txt";
	p->myFile = theDefault;
	//int i =0;
	//while(p->myFile[i] = theDefault[i]){
	//  i++;
	//}
	std::cout << "The file name has been changed to: \"" << p->myFile << "\"\n";
	//{ p->myFile = "../ACOorbits.txt";// initalize a default file
      }*/
  }
    else if(v.varName == "filename") {
      p->myFile = v.varString;
      std::cout << "The file name has been changed to: \"" << p->myFile << "\"\n";
    }
    else if(v.varName == "timeStep") {
      p->timeStep = v.varVal*kUnits::sec;
      std::cout << "timeStep is " << p->fileMode/kUnits::s << " seconds"  << std::endl;
    }
    else if(v.varName == "excitedOnly") {
      p->excitedOnly = bool(v.varVal);
      //p->solarNeutrinoRate = User::neutrinoRate(0,20)*p->galliumAtomNumber*0.399; // update rate
      // rate is updated in the main program because it's clunky here
      std::cout << "excitedOnly is " << p->excitedOnly <<  std::endl;
    }
    else if(v.varName == "detectorRadius") {
      p->detectorRadius = v.varVal*kUnits::meter;
      std::cout << "detectorRadius is " << p->detectorRadius/kUnits::meter << " meters"  << std::endl;
    }
    else if(v.varName == "kgGallium") {
      p->kgGallium = v.varVal*kUnits::kg;
      p->galliumAtomNumber = p->kgGallium * 8.637e24/kUnits::kg;// calculate the derived atom number   
      std::cout << "kgGallium is " << p->kgGallium/kUnits::kg << " kg, and galliumAtomNumber is "  
		<< p->galliumAtomNumber << std::endl;
      // recalculate the new rate after updating the kg of gallium
      //p->solarNeutrinoRate = User::neutrinoRate(0,20)*p->galliumAtomNumber*0.399;
      // rate is updated in the main program because it's clunky here
    }
    else if(v.varName == "closest") {
      p->closest = v.varVal*kUnits::solarRadii;
      std::cout << "closest is " << p->closest/kUnits::solarRadii << " solar radii"  << std::endl;
    }
    else if(v.varName == "furthest") {
      p->furthest = v.varVal*kUnits::solarRadii;
      std::cout << "furthest is " << p->furthest/kUnits::solarRadii << " solar radii"  << std::endl;
    }
    else if(v.varName == "nLoops") {
      p->nLoops = v.varVal;
      std::cout << "nLoops is " << p->nLoops << std::endl;
      if( p->timeLimit > 0){
      std::cout << "timeLimit is positive. Setting nLoops to 4294967295" << std::endl;
      p->nLoops = 4294967295;
      }
    }
    else if(v.varName == "timeLimit") {
      p->timeLimit = v.varVal * kUnits::year;
      std::cout << "timeLimit (in days) is " << p->timeLimit/kUnits::day << std::endl;
      if( p->timeLimit > 0){
      std::cout << "timeLimit is positive. Setting nLoops to 4294967295" << std::endl;
      p->nLoops = 4294967295;
      }
    }
    else{
      std::cout << "This variable wasn't one I know.\n";
    }
  }


  char* getCmdOption(char ** begin, char ** end, const std::string & option)
  {
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
      {
        return *itr;
      }
    return 0;
  }

  bool cmdOptionExists(char** begin, char** end, const std::string& option)
  {
    return std::find(begin, end, option) != end;
  }

  bool isComment(std::string theString)
  {
    return theString.at(0) == '#';
  }

  void variablePusher(std::vector<User::Variable> &holder, std::string file){
    // These will hold the contents and each line
    std::ifstream fileContents;
    std::string line;

    // open the file
    fileContents.open(file);

    while(getline(fileContents, line)){
      if(isComment(line));// do nothing for comments
      else{
	// take line and make it into an sstream (pieces separated by ' '
	std::stringstream words(line);
	// set up holding variables
	std::string paramName;
	std::string paramString;
	User::Variable toPush;
	// fill holding variables
	words >> paramName >> paramString;
	double paramValue;
	toPush.varName = paramName;
	toPush.varString = paramString;
	try {
	  paramValue = std::stod(paramString);
	  toPush.varVal = paramValue;
	  
	} catch (const std::invalid_argument& e) {
	  //std::cerr << "Invalid argument for stod: " << paramString << std::endl;
	} catch (const std::out_of_range& e) {
	  //std::cerr << "Argument out of range for stod: " << paramString << std::endl;
	}
	
	/*
	if(paramValue = stod(paramString)){
	  toPush.varVal = paramValue;
	}
	*/
	// push onto vector
	holder.push_back(toPush);
	std::cout << "I have pushed " << paramString << " into " << paramName 
		  << " and the holding vector is " << holder.size() 
		  << " elements long\n";
      }
    }
  
  
  
  }

  void printHelp(){
    printedHelp = true;
    std::cout << "\n\nUsage: ./nuSolPerformance <options>\n\n  Options are:\n\n    -h : Prints the help screen\n    -b : Run in batch mode. Follow with the batch file. E.g. ./nuSolPerformance -b batch.mac\n    -f : DEPRICATED provide a filepath for the flight input (put in batch file)\n\n\n";

  }
  
  int argumentInterpreter(int argc, char * argv[], User::Par *par)
  {
    parInitializer(par); // initialize base values for arguments
    if(cmdOptionExists(argv, argv+argc, "-h"))
      {
	printHelp();
      }

    char * filename = getCmdOption(argv, argv + argc, "-b");

    if (filename)
      {
	std::vector<User::Variable> V;
	variablePusher(V,filename);
	User::Variable v;
	std::cout << "I have " << V.size() << " variables to change.\n";
	while(!V.empty()){
	  v = V.back();
	  parChanger(par,v);
	  V.pop_back();
	}
        // Do interesting things
        // ...
      }

    char * filename1 = getCmdOption(argv, argv + argc, "-f");

    if (filename1)
      {
	std::string s(filename1);
	char* theString = new char[s.length()+1];
	for(int i=0; i<s.length(); i++){
	  par->myFile[i] = s[i];
	  //std::cout << "I just wrote '" << par->myFile[i] << "' to the filename from the original value of '" << s[i] << "'\n";
	}
	std::cout << "The position file to open is " << par->myFile << '\n';
      }

    return 0;
  }









}




// First we read in the file we're sent from the main function
