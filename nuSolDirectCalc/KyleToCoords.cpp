#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

int main() {
  std::ifstream inFile("KyleOrbit.csv");
  std::ofstream outFile("KyleOrbitParsed.csv");

  if (!inFile.is_open()) {
    std::cerr << "Failed to open input file." << std::endl;
    return 1;
  }

  if (!outFile.is_open()) {
    std::cerr << "Failed to open output file." << std::endl;
    return 1;
  }

  std::string line;
  double lastValue = 0.0;
  bool firstLine = true;

  while (getline(inFile, line)) {
    std::stringstream lineStream(line);
    std::string cell;
    double currentValue = 0.0;
    int columnIndex = 0;

    while (getline(lineStream, cell, ',')) {
      if (columnIndex == 0) {
	currentValue = std::stod(cell);
	if (firstLine || (currentValue - lastValue > 0)) {
	  outFile << line << '\n';
	  lastValue = currentValue;
	  firstLine = false;
	  break; // No need to process the rest of the line
	}
      }
      columnIndex++;
    }
  }

  inFile.close();
  outFile.close();

  return 0;
}
