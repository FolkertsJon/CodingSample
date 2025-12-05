#include <iostream>
#include <fstream>

int main() {
    // Open the file in write mode which overwrites any existing content
    std::ofstream outFile("lineOrbit.csv", std::ofstream::out);

    // Check if the file is open
    if (!outFile.is_open()) {
        std::cerr << "Failed to open file." << std::endl;
        return 1;
    }

    double au = 149597870700;// meters
    double RSol = 695700000;// meters

    // Iterate from 0 to 99 with 0.01 steps (Moves from 100 to 1 solar radii)
    for (double t = 0; t <= 99; t += 0.01) {
      outFile << t << "," << (100 - t)*RSol/au << "\n";
    }

    // Close the file
    outFile.close();

    return 0;
}
