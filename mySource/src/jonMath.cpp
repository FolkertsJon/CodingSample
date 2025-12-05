// This program is designed to allow for several oscillation models
// to be implemented in the main file. These models will output
// a multiplicative factor to be used in the main program which
// is generally a function of radius.


#include<jonMath.hh>

namespace User{
  
  double trapezoidIntegral(double xMin, double yMin, double xMax, double yMax, double intMin, double intMax){
    
    // restrict to the range of interest
    double theMin = std::max(intMin,xMin);
    double theMax = std::min(intMax,xMax);
    double integral = (theMax - theMin) * (yMax + yMin)/2;
    return integral;
  }

  double trapezoidIntegral(double xMin, double yMin, double xMax, double yMax){
    return trapezoidIntegral(xMin, yMin, xMax, yMax, xMin, xMax);
  }

  double trapezoidIntegral(std::vector<double> x, std::vector<double> y){
    double integral = 0;
    for( int i = 1; i < x.size(); i++){// skip first for numerical integration
      double toAdd = trapezoidIntegral(x[i-1],y[i-1],x[i],y[i]);
      integral += toAdd;
    }
  return integral;
  }

  
  double trapezoidIntegral(std::vector<double> x, std::vector<double> y, double xMax, double yMax){
    double integral = 0;
    for( int i = 1; i < x.size(); i++){// skip first for numerical integration
      double toAdd = trapezoidIntegral(x[i-1],y[i-1],x[i],y[i], xMax, yMax);
      integral += toAdd;
    }
  return integral;
  }


  std::pair<long double, long double> wilsonScoreInterval(size_t count, size_t total, long double z){
    if (total == 0) return {0.0, 0.0}; // To avoid division by zero
    long double p = static_cast<long double>(count) / total;
    long double denominator = 1.0 + (z * z) / total;
    long double margin = z * std::sqrt((p * (1 - p) / total) + (z * z) / (4.0 * total * total));
    long double lower = (p + (z * z) / (2.0 * total) - margin) / denominator;
    long double upper = (p + (z * z) / (2.0 * total) + margin) / denominator;
    return {lower, upper};
  }


  double solveKepler(double M, double e) {
    double E = M, E1;
    double tolerance = 1e-6; // tolerance for Newton's method convergence
    int maxIter = 100, i = 0;

    do {
      double f = E - e * sin(E) - M;
      double df = 1 - e * cos(E);
      E1 = E - f / df;
      if (fabs(E1 - E) < tolerance && i > 9) // do at least 10 iterations
	break;
      E = E1;
      i++;
    } while (i < maxIter);

    if (i == maxIter) {
      std::cout << "Maximum iterations reached, returning answer.\n";
    }
    return E1;
  }

  double computePositionForKepler(double a, double e, double E) {
    // results are identical to the ChatGPT version to 15 digits in one test orbit.
    double b = a * sqrt(1 - e*e); // calculate semi-minor axis
    double x = a * (cos(E) - e); // x from Kepler
    double y = b * sin(E); // y from Kepler
    return sqrt(x*x + y*y); // Calculate radius and return
  }

  double computePositionForKeplerChatGPT(double a, double e, double E) {
    return a * (1 - e * cos(E)); // Calculate radius and return
  }



  double betheBlochLoss(double z, double a, double beta, double density){
    // Constants from kUnits, explicitly qualified
    const double me = kUnits::me;        // Electron mass
    const double e = kUnits::elementaryCharge;     // Elementary charge
    const double eps0 = kUnits::vacuumPermittivity;// Vacuum permittivity
    const double c = kUnits::speedOfLight;         // Speed of light

    std::cout << "c = " << c*kUnits::s/kUnits::m << " m/s\n";
    std::cout << "eps0 = " << eps0/kUnits::coulomb/kUnits::coulomb/kUnits::s/kUnits::s*kUnits::kg*kUnits::m*kUnits::m*kUnits::m << "\n";
    
    double I = 10*kUnits::eV*z;
    // Number density of electrons
    // n = Z/A * N_A * density
    double n = z * kUnits::AvogadroNumber * density / a /kUnits::gram;
    std::cout << "n = " << n*kUnits::m*kUnits::m*kUnits::m << "/m^3\n";

    // The constant factor K, with units squared to match the usage
    // double K = 4 * M_PI / (me * c *c) * n/(beta*beta)*(e*e/(4*M_PI*eps0)*(4*M_PI*eps0));
    double term1 =  4 * M_PI / (me * c *c);
    std::cout << "term1 = " << term1*kUnits::kg/(kUnits::m/kUnits::s)/(kUnits::m/kUnits::s) << " m^2 /kg /s^2\n";
    double term2 =  n*z*z/(beta*beta);
    std::cout << "term2 = " << term2 * kUnits::m * kUnits::m * kUnits::m << " /m^3\n";
    double term3 = (e*e/(4*M_PI*eps0))*(e*e/(4*M_PI*eps0));
    double term3Units = term3 / kUnits::kg / kUnits::kg / kUnits::m / kUnits::m / kUnits::m / kUnits::m / kUnits::m / kUnits::m * kUnits::s * kUnits::s * kUnits::s * kUnits::s;
    std::cout << "term3  = " << term3Units << " kg^2 m^6 s^-4\n";
    
    // Bethe formula
    double gamma = 1.0 / sqrt(1.0 - beta * beta);
    double Emax = 2.0 * me * std::pow(c, 2) * beta * beta * gamma * gamma;
    double lnTerm = std::log(Emax / I) - beta*beta;
    std::cout << "lnTerm = " << lnTerm << "\n";

    // Energy loss per unit length
    double dEdx = term1 * term2 * term3 * lnTerm;

    return dEdx;
  }
  
}

