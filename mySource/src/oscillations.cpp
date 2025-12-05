// This program is designed to allow for several oscillation models
// to be implemented in the main file. These models will output
// a multiplicative factor to be used in the main program which
// is generally a function of radius.


#include<oscillations.hh>

namespace User
{

  std::vector<User::energyProbPair> readData(const std::string& filename) {
    std::ifstream inFile(filename);
    if (!inFile) {
      std::cerr << "Unable to open file " << filename << "\n";
      exit(1);
    }

    std::vector<User::energyProbPair> data;
    std::string line;
    while (std::getline(inFile, line)) {
      std::stringstream ss(line);
      double e, p;
      if (ss >> e >> p) {
	data.push_back(std::make_pair(e, p));
      }
    }
    inFile.close();
    
    // Ensure the data is sorted by energy for interpolation
    std::sort(data.begin(), data.end());

    return data;
  }
  
  // thic returns the probabiltiy of a neutrino remaining electron-type
  // at the earth when coming from the sun.
  double nu_e_Earth(const std::vector<User::energyProbPair>& data, double energy) {
    // Extrapolation below the range
    if (energy < data.front().first) {
        double slope = (data[1].second - data[0].second) / (data[1].first - data[0].first);
        return data[0].second + slope * (energy - data[0].first);
    }

    // Extrapolation above the range
    if (energy > data.back().first) {
        double slope = (data.back().second - data[data.size() - 2].second) / (data.back().first - data[data.size() - 2].first);
        return data.back().second + slope * (energy - data.back().first);
    }

    // Regular case: within the range
    for (size_t i = 0; i < data.size() - 1; ++i) {
        if (data[i].first <= energy && energy <= data[i+1].first) {
            // Linear interpolation
            double slope = (data[i+1].second - data[i].second) / (data[i+1].first - data[i].first);
            return data[i].second + slope * (energy - data[i].first);
        }
    }

    // Unreachable code, since we've handled all cases
    return -1.0;
}




}
