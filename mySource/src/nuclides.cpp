

// This source file is designe for any custom nuclide info I need to
// parse

#include "nuclides.hh"

namespace User{
  // Function to map Z numbers to element symbols
  std::string getElementSymbol(int Z) {
    static const std::unordered_map<int, std::string> symbolMap = {
      {0,"freeNeutron"}, {1, "H"}, {2, "He"}, {3, "Li"}, {4, "Be"}, {5, "B"},
      {6, "C"}, {7, "N"}, {8, "O"}, {9, "F"}, {10, "Ne"},
      {11, "Na"}, {12, "Mg"}, {13, "Al"}, {14, "Si"}, {15, "P"},
      {16, "S"}, {17, "Cl"}, {18, "Ar"}, {19, "K"}, {20, "Ca"},
      {21, "Sc"}, {22, "Ti"}, {23, "V"}, {24, "Cr"}, {25, "Mn"},
      {26, "Fe"}, {27, "Co"}, {28, "Ni"}, {29, "Cu"}, {30, "Zn"},
      {31, "Ga"}, {32, "Ge"}, {33, "As"}, {34, "Se"}, {35, "Br"},
      {36, "Kr"}, {37, "Rb"}, {38, "Sr"}, {39, "Y"}, {40, "Zr"},
      {41, "Nb"}, {42, "Mo"}, {43, "Tc"}, {44, "Ru"}, {45, "Rh"},
      {46, "Pd"}, {47, "Ag"}, {48, "Cd"}, {49, "In"}, {50, "Sn"},
      {51, "Sb"}, {52, "Te"}, {53, "I"}, {54, "Xe"}, {55, "Cs"},
      {56, "Ba"}, {57, "La"}, {58, "Ce"}, {59, "Pr"}, {60, "Nd"},
      {61, "Pm"}, {62, "Sm"}, {63, "Eu"}, {64, "Gd"}, {65, "Tb"},
      {66, "Dy"}, {67, "Ho"}, {68, "Er"}, {69, "Tm"}, {70, "Yb"},
      {71, "Lu"}, {72, "Hf"}, {73, "Ta"}, {74, "W"}, {75, "Re"},
      {76, "Os"}, {77, "Ir"}, {78, "Pt"}, {79, "Au"}, {80, "Hg"},
      {81, "Tl"}, {82, "Pb"}, {83, "Bi"}, {84, "Po"}, {85, "At"},
      {86, "Rn"}, {87, "Fr"}, {88, "Ra"}, {89, "Ac"}, {90, "Th"},
      {91, "Pa"}, {92, "U"}, {93, "Np"}, {94, "Pu"}, {95, "Am"},
      {96, "Cm"}, {97, "Bk"}, {98, "Cf"}, {99, "Es"}, {100, "Fm"},
      {101, "Md"}, {102, "No"}, {103, "Lr"}, {104, "Rf"}, {105, "Db"},
      {106, "Sg"}, {107, "Bh"}, {108, "Hs"}, {109, "Mt"}, {110, "Ds"},
      {111, "Rg"}, {112, "Cn"}, {113, "Nh"}, {114, "Fl"}, {115, "Mc"},
      {116, "Lv"}, {117, "Ts"}, {118, "Og"}
    };

    auto it = symbolMap.find(Z);
    if (it != symbolMap.end()) {
      return it->second;
    } else {
      std::cerr << "Element Z not found: " << Z << std::endl;
      //return -1; // Error code for not found
    }
  }

  
  int getAtomicNumber(const std::string& symbol) {
    static const std::unordered_map<std::string, int> atomicNumberMap = {
      {"freeNeutron", 0}, {"H", 1}, {"He", 2}, {"Li", 3}, {"Be", 4}, {"B", 5},
      {"C", 6}, {"N", 7}, {"O", 8}, {"F", 9}, {"Ne", 10},
      {"Na", 11}, {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"P", 15},
      {"S", 16}, {"Cl", 17}, {"Ar", 18}, {"K", 19}, {"Ca", 20},
      {"Sc", 21}, {"Ti", 22}, {"V", 23}, {"Cr", 24}, {"Mn", 25},
      {"Fe", 26}, {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30},
      {"Ga", 31}, {"Ge", 32}, {"As", 33}, {"Se", 34}, {"Br", 35},
      {"Kr", 36}, {"Rb", 37}, {"Sr", 38}, {"Y", 39}, {"Zr", 40},
      {"Nb", 41}, {"Mo", 42}, {"Tc", 43}, {"Ru", 44}, {"Rh", 45},
      {"Pd", 46}, {"Ag", 47}, {"Cd", 48}, {"In", 49}, {"Sn", 50},
      {"Sb", 51}, {"Te", 52}, {"I", 53}, {"Xe", 54}, {"Cs", 55},
      {"Ba", 56}, {"La", 57}, {"Ce", 58}, {"Pr", 59}, {"Nd", 60},
      {"Pm", 61}, {"Sm", 62}, {"Eu", 63}, {"Gd", 64}, {"Tb", 65},
      {"Dy", 66}, {"Ho", 67}, {"Er", 68}, {"Tm", 69}, {"Yb", 70},
      {"Lu", 71}, {"Hf", 72}, {"Ta", 73}, {"W", 74}, {"Re", 75},
      {"Os", 76}, {"Ir", 77}, {"Pt", 78}, {"Au", 79}, {"Hg", 80},
      {"Tl", 81}, {"Pb", 82}, {"Bi", 83}, {"Po", 84}, {"At", 85},
      {"Rn", 86}, {"Fr", 87}, {"Ra", 88}, {"Ac", 89}, {"Th", 90},
      {"Pa", 91}, {"U", 92}, {"Np", 93}, {"Pu", 94}, {"Am", 95},
      {"Cm", 96}, {"Bk", 97}, {"Cf", 98}, {"Es", 99}, {"Fm", 100},
      {"Md", 101}, {"No", 102}, {"Lr", 103}, {"Rf", 104}, {"Db", 105},
      {"Sg", 106}, {"Bh", 107}, {"Hs", 108}, {"Mt", 109}, {"Ds", 110},
      {"Rg", 111}, {"Cn", 112}, {"Nh", 113}, {"Fl", 114}, {"Mc", 115},
      {"Lv", 116}, {"Ts", 117}, {"Og", 118}
    };

    auto it = atomicNumberMap.find(symbol);
    if (it != atomicNumberMap.end()) {
      return it->second;
    } else {
      std::cerr << "Element symbol not found: " << symbol << std::endl;
      return -1; // Error code for not found
    }
  }



  int getAtomicNumber(int Z){
    return getAtomicNumber(getElementSymbol(Z));
  }
  

  double getAtomicMass(const std::string& symbol) {
    static const std::unordered_map<std::string, double> atomicNumberMap = {
      {"freeNeutron", 0}, {"H", 1.0080}, {"He", 4.00260}, {"Li", 7.0},
      {"Be", 9.012183}, {"B", 10.81}, {"C", 12.011}, {"N", 14.007},
      {"O", 15.999}, {"F", 18.99840316}, {"Ne", 20.180}, {"Na", 22.9897693 },
      {"Mg", 24.305}, {"Al", 26.981538}, {"Si", 28.085}, {"P", 30.97376200},
      {"S", 32.07}, {"Cl", 35.45}, {"Ar", 39.9}, {"K", 39.0983}, {"Ca", 40.08},
      {"Sc", 44.95591}, {"Ti", 47.867}, {"V", 50.9415}, {"Cr", 51.996},
      {"Mn", 54.93804}, {"Fe", 55.84}, {"Co", 58.93319}, {"Ni", 58.693},
      {"Cu", 63.55}, {"Zn", 65.4}, {"Ga", 69.723}, {"Ge", 72.63},
      {"As", 74.92159}, {"Se", 78.97}, {"Br", 79.90}, {"Kr", 83.80},
      {"Rb", 85.468}, {"Sr", 87.62}, {"Y", 88.90584}, {"Zr", 91.22},
      {"Nb", 92.90637}, {"Mo", 95.95}, {"Tc", 96.90636}, {"Ru", 101.1},
      {"Rh", 102.9055}, {"Pd", 106.42}, {"Ag", 107.86}, {"Cd", 112.41},
      {"In", 114.818}, {"Sn", 118.71}, {"Sb", 121.760}, {"Te", 127.6},
      {"I", 126.9045}, {"Xe", 131.29}, {"Cs", 132.9054520}, {"Ba", 137.33},
      {"La", 138.9055}, {"Ce", 140.116}, {"Pr", 140.90766}, {"Nd", 144.24},
      {"Pm", 144.91276}, {"Sm", 150.4}, {"Eu", 151.964}, {"Gd", 157.2},
      {"Tb", 158.92535}, {"Dy", 162.500}, {"Ho", 164.93033}, {"Er", 167.26},
      {"Tm", 168.93422}, {"Yb", 173.05}, {"Lu", 174.9668}, {"Hf", 178.49},
      {"Ta", 180.9479}, {"W", 183.84}, {"Re", 186.207}, {"Os", 190.2},
      {"Ir", 192.22}, {"Pt", 195.08}, {"Au", 196.96657}, {"Hg", 200.59},
      {"Tl", 204.383 }, {"Pb", 207}, {"Bi", 208.98040}, {"Po", 208.98243},
      {"At", 209.98715}, {"Rn", 222.01758}, {"Fr", 223.01973},
      {"Ra", 226.02541}, {"Ac", 227.02775}, {"Th", 232.038}, {"Pa", 231.03588},
      {"U", 238.0289}, {"Np", 237.048172}, {"Pu", 244.06420},
      {"Am", 243.061380}, {"Cm", 247.07035}, {"Bk", 247.07031},
      {"Cf", 251.07959}, {"Es", 252.0830}, {"Fm", 257.09511}, {"Md", 258.09843},
      {"No", 259.10100}, {"Lr", 266.120}, {"Rf", 267.122}, {"Db", 268.126},
      {"Sg", 269.128}, {"Bh", 270.133}, {"Hs", 269.1336}, {"Mt", 277.154},
      {"Ds", 282.166}, {"Rg", 282.169}, {"Cn", 286.179}, {"Nh", 286.182},
      {"Fl", 290.192}, {"Mc", 290.196 }, {"Lv", 293.205}, {"Ts", 294.211},
      {"Og", 295.216}
    };

    auto it = atomicNumberMap.find(symbol);
    if (it != atomicNumberMap.end()) {
      return it->second;
    } else {
      std::cerr << "Element symbol not found: " << symbol << std::endl;
      return -1; // Error code for not found
    }
  }
  
  double getAtomicMass(int Z){
    return getAtomicMass(getElementSymbol(Z));
  }
  
  std::pair<int, int> ZAFromPDGCode(const std::string& pdgId){
    if (pdgId.length() != 10) {
      throw std::invalid_argument("PDG ID must be a 10-digit string.");
    }

    // Extracting Z and A
    std::string zStr = pdgId.substr(3, 3);  // Z is located at positions 3, 4, 5
    std::string aStr = pdgId.substr(6, 3);  // A is located at positions 6, 7, 8

    // Converting string to integer
    int Z = std::stoi(zStr);
    int A = std::stoi(aStr);

    return std::make_pair(Z, A);
  }

  
}// End User Namespace
