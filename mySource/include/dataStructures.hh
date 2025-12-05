// This holds the custom data structures I use for my functions
// Jonathan Folkerts June 2021
#ifndef dataStructures_H
#define dataStructures_H

#include<string>
#include<vector>
#include<complex>
#include<Eigen/Dense>

#include<kUnits.hh>

namespace User
{
  struct Par{
    // Control Variables!
    bool fileMode;
    double timeStep;// time step for iterating
    //char myFile[50]; // String that holds r vs time
    std::string myFile;
    double timeLimit; // soft cap on execution time
    double detectorRadius; // approximate value
    double cosmicAcceptance; // 5 nines
    double solarAcceptance;// four nines better than cosmic
    double kgGallium;
    double galliumAtomNumber;//
    double closest; // Closest approach for elliptical mode
    double furthest; // Furthest approach for elliptical mode
    uint nLoops; // number of iterations to run
    uint nSteps; // Used by file mode
    double solarNeutrinoRate; // per atom rate * nAtoms
    double cosBack; // cosmic background rate
    double inclination; // currently unused
    bool excitedOnly; // boolean in readout
  };

  struct HistoDat{
    // Variables to put into treeeeeees
    double neutrinoRadius = 0;
    double sub35NeutrinoRadius = 0;
    //double neutrinoOnlyRadius;
    //double sub35NeutrinoOnlyRadius;
    double cosmicBackRadius = 0;
    double sub35CosmicBackRadius = 0;
    double solarBackRadius = 0;
    double sub35SolarBackRadius = 0;
    double radius = 0;
  };


  struct csvOutputs {
    std::vector<double> timeOfRise;
    std::vector<double> timeOfMax;
    std::vector<double> digitalTimes;
    std::vector<int> maxIdx;
    std::vector<double> voltagePeaks;
    std::vector<double> voltageIntegrations;
    std::vector<double> channelThresholds;
    std::vector<int> properSign;
    std::vector<bool> risen;

    // Explicitly delete the default constructor
    csvOutputs() = delete;

    // Constructor that sets 'channelThresholds' to the provided vector
    // and initializes time vectors with NaN values
    csvOutputs(const std::vector<double>& vec) 
      : timeOfRise(vec.size(), std::nan("")), 
	timeOfMax(vec.size(), std::nan("")), 
	digitalTimes(16, std::nan("")), 
	maxIdx(vec.size()), 
	voltagePeaks(vec.size(), 0), 
	voltageIntegrations(vec.size(), 0), 
	channelThresholds(vec), 
	properSign(vec.size(), 0), 
	risen(vec.size(), false)  
    {}
  };

  
  struct Variable{
    std::string varName;
    double varVal;
    std::string varString;
  };// This holds a variable name and a variable value for changing parameters



  
  struct neutrinoData{
    // Probability holder for neutrino flavors (0=electron, 1=mu, 2=tau)
    double eProb;
    double muProb;
    double tauProb;
    double eDev;
    double muDev;
    double tauDev;
    
    // neutrino Ket
    Eigen::Vector3cd nuKet;
    
    // neutrino energy (approximate; it's for averaging
    double neutrinoEnergy;
    
    // neutrino current position in solar radii
    double X;
    double Y;
    double Z;
    double neutrinoPosition; // r in solar radii

    // distance moved in last movement step
    double justMoved;

    // distance of next shell in sun model and step size
    double rNext;
    double rStep;
    
    // neutrino starting position in solar radii
    double neutrinoStartR;
    double neutrinoStartPhi;
    double neutrinoStartTheta;

  };

  struct neutrinoTarget {
    std::string name;        // Target's name
    int PDGcode = 0;         // Particle Data Group code
    double amuPerAtom = 0.0; // Atomic mass unit per atom
    int atomicNumber = 0;    // Atomic number
    int massNumber = 0;      // Mass number
    double nPerKilogram = 0.0; // number of targets per kilogram
    double nPerMass = 0.0; // number of targets permass
    double massPerAtom = 0.0;  // mass per atom
    double targetMassFraction = 0.0;  // mass fraction of desired target

    // Default constructor
    neutrinoTarget() 
      : name(""), PDGcode(0), amuPerAtom(0.0), atomicNumber(0), massNumber(0),
	nPerKilogram(0.0), massPerAtom(0.0), targetMassFraction(0.0) {}
  };
  
  // First is energy, second is Probability of remaing electron-type neutrino
  typedef std::pair<double, double> energyProbPair;

  
  struct oscilloscopeData {
    std::string oscilloscopeModel; 
    std::string horizontalUnit; 
    double horizontalScale;
    double horizontalDigitalScale;
    std::vector<size_t> recordLengths; // how long each channel is; constant for tek Oscope
    std::vector<size_t> digitalRecordLengths; // how long each channel is; constant for tek Oscope
    std::vector<double> probeAttenuation; // attenuation for each channel
    std::vector<std::string> verticalUnits;
    std::vector<std::string> digitalVerticalUnits;
    std::vector<double> verticalScale;
    std::vector<double> verticalMinorTic;
    std::vector<double> verticalPosition;
    std::vector< std::vector<double> > channelData;
    std::vector<double> timeData;
    std::vector< std::vector<double> > digitalData;
    std::vector<double> digitalTimeData;
    std::vector<std::string> channelLabels;
    size_t numberOfChannels = 0;
    size_t numberOfDigitalChannels = 0;

    friend std::ostream& operator<<(std::ostream& os, const oscilloscopeData& data);
  };

  inline std::ostream& operator<<(std::ostream& os, const oscilloscopeData& data) {
    os << "oscilloscopeModel: " << data.oscilloscopeModel << "\n";
    os << "horizontalUnit: " << data.horizontalUnit << "\n";
    os << "horizontalScale: " << data.horizontalScale << "\n";
        
    os << "recordLengths: ";
    for (size_t length : data.recordLengths) {
      os << length << " ";
    }
    os << "\n";

    os << "probeAttenuation: ";
    for (double attenuation : data.probeAttenuation) {
      os << attenuation << " ";
    }
    os << "\n";

    os << "verticalUnits: ";
    for (const std::string& unit : data.verticalUnits) {
      os << unit << " ";
    }
    os << "\n";

    os << "verticalScale: ";
    for (double scale : data.verticalScale) {
      os << scale << " ";
    }
    os << "\n";

    os << "verticalMinorTic: ";
    for (double tic : data.verticalMinorTic) {
      os << tic << " ";
    }
    os << "\n";

    os << "verticalPosition: ";
    for (double position : data.verticalPosition) {
      os << position << " ";
    }
    os << "\n";

    os << "numberOfChannels: " << data.numberOfChannels << "\n";
    return os;
  }
}

#endif
