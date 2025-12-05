/*
  This is a program to test using Eigen and to get the neutrino oscillation
  calcluations underway. The wave propagates as e^-iHt/hbar, so we need to
  calculate 1+H+H^2/2!+... until there is less than x% change

  This file now finds the probability of being in all three of the neutrino states
  
*/

#include "neutrinoOscillation.hh"

namespace User
{
  
  Eigen::Matrix3cd hamiltonianStep(User::neutrinoData *theData){
    double Y1;
    if(theData -> Y < 0){
      Y1 = -(theData -> Y);// move to y=0 plane
    }
    else{
      Y1 = -(theData -> Y) + sqrt( (theData -> Y)*(theData -> Y)
		      +(theData -> rNext)*(theData -> rNext)
		      -(theData -> X)*(theData -> X)
		      -(theData -> Z)*(theData -> Z) 
		     );// move to next ring
      //std::cout << "The distance to move is: " << Y1/kUnits::solarRadii << " solar radii.\n"
    }
    
    double rOld = sqrt( (theData -> Y)*(theData -> Y)
			+(theData -> X)*(theData -> X)
			+(theData -> Z)*(theData -> Z));// find old R
    double rNew = sqrt( ((theData -> Y) + Y1)*((theData -> Y) + Y1)
			+(theData -> X)*(theData -> X)
			+(theData -> Z)*(theData -> Z));//find new R

    Eigen::Matrix3cd theReturn = sunSectionHamiltonian ( std::min(rOld,rNew), std::max(rOld,rNew) , theData ->neutrinoEnergy);// calculate the hamiltonian to return

    theData -> justMoved = Y1;// distance we just moved
    theData -> Y += Y1 ;// move to the next sphere (+1 meter for safety)
    theData -> rNext += theData -> rStep;
    
    return theReturn;
  }

  Eigen::Matrix3cd hamiltonianStepIO(User::neutrinoData *theData){
    double Y1;
    if(theData -> Y < 0){
      Y1 = -(theData -> Y);// move to y=0 plane
    }
    else{
      Y1 = -(theData -> Y) + sqrt( (theData -> Y)*(theData -> Y)
		      +(theData -> rNext)*(theData -> rNext)
		      -(theData -> X)*(theData -> X)
		      -(theData -> Z)*(theData -> Z) 
		     );// move to next ring
      //std::cout << "The distance to move is: " << Y1/kUnits::solarRadii << " solar radii.\n"
    }
    
    double rOld = sqrt( (theData -> Y)*(theData -> Y)
			+(theData -> X)*(theData -> X)
			+(theData -> Z)*(theData -> Z));// find old R
    double rNew = sqrt( ((theData -> Y) + Y1)*((theData -> Y) + Y1)
			+(theData -> X)*(theData -> X)
			+(theData -> Z)*(theData -> Z));//find new R

    Eigen::Matrix3cd theReturn = sunSectionHamiltonianIO ( std::min(rOld,rNew), std::max(rOld,rNew) , theData ->neutrinoEnergy);// calculate the hamiltonian to return

    theData -> justMoved = Y1;// distance we just moved
    theData -> Y += Y1 ;// move to the next sphere (+1 meter for safety)
    theData -> rNext += theData -> rStep;
    
    return theReturn;
  }

  Eigen::Matrix3cd U(){//3-d complex double matrix
    
    std::complex<double> U11 = cos(kUnits::theta12)*cos(kUnits::theta13);
    std::complex<double> U12 = sin(kUnits::theta12)*cos(kUnits::theta13);
    std::complex<double> U13 = sin(kUnits::theta13)*std::exp(-kUnits::i * kUnits::deltacp);
    std::complex<double> U21 =
      -sin(kUnits::theta12)*cos(kUnits::theta23) - cos(kUnits::theta12)*sin(kUnits::theta23)*sin(kUnits::theta13)*std::exp(kUnits::i * kUnits::deltacp);
    std::complex<double> U22 =
      cos(kUnits::theta12)*cos(kUnits::theta23) - sin(kUnits::theta12)*sin(kUnits::theta23)*sin(kUnits::theta13)*std::exp(kUnits::i * kUnits::deltacp);
    std::complex<double> U23 = sin(kUnits::theta23)*cos(kUnits::theta13);
    std::complex<double> U31 =
      sin(kUnits::theta12)*sin(kUnits::theta23)-cos(kUnits::theta12)*cos(kUnits::theta23)*sin(kUnits::theta13)*std::exp(kUnits::i * kUnits::deltacp);
    std::complex<double> U32 =
      -cos(kUnits::theta12)*sin(kUnits::theta23)-sin(kUnits::theta12)*cos(kUnits::theta23)*sin(kUnits::theta13)*std::exp(kUnits::i * kUnits::deltacp);
    std::complex<double> U33 = cos(kUnits::theta23)*cos(kUnits::theta13);
  
    Eigen::Matrix3cd theReturn;
    theReturn << U11, U12, U13,
      U21, U22, U23,
      U31, U32, U33;
    return theReturn;
  }

  Eigen::Matrix3cd UIO(){//3-d complex double matrix
    
    std::complex<double> U11 = cos(kUnits::theta12IO)*cos(kUnits::theta13IO);
    std::complex<double> U12 = sin(kUnits::theta12IO)*cos(kUnits::theta13IO);
    std::complex<double> U13 = sin(kUnits::theta13IO)*std::exp(-kUnits::i * kUnits::deltacpIO);
    std::complex<double> U21 =
      -sin(kUnits::theta12IO)*cos(kUnits::theta23IO) - cos(kUnits::theta12IO)*sin(kUnits::theta23IO)*sin(kUnits::theta13IO)*std::exp(kUnits::i * kUnits::deltacpIO);
    std::complex<double> U22 =
      cos(kUnits::theta12IO)*cos(kUnits::theta23IO) - sin(kUnits::theta12IO)*sin(kUnits::theta23IO)*sin(kUnits::theta13IO)*std::exp(kUnits::i * kUnits::deltacpIO);
    std::complex<double> U23 = sin(kUnits::theta23IO)*cos(kUnits::theta13IO);
    std::complex<double> U31 =
      sin(kUnits::theta12IO)*sin(kUnits::theta23IO)-cos(kUnits::theta12IO)*cos(kUnits::theta23IO)*sin(kUnits::theta13IO)*std::exp(kUnits::i * kUnits::deltacpIO);
    std::complex<double> U32 =
      -cos(kUnits::theta12IO)*sin(kUnits::theta23IO)-sin(kUnits::theta12IO)*cos(kUnits::theta23IO)*sin(kUnits::theta13IO)*std::exp(kUnits::i * kUnits::deltacpIO);
    std::complex<double> U33 = cos(kUnits::theta23IO)*cos(kUnits::theta13IO);
  
    Eigen::Matrix3cd theReturn;
    theReturn << U11, U12, U13,
      U21, U22, U23,
      U31, U32, U33;
    return theReturn;
  }

  Eigen::Matrix3cd k(){// k matrix for mixing; normal ordering
    Eigen::Matrix3cd theReturn;
    theReturn << 0, 0, 0,
      0, kUnits::deltam21, 0,
      0, 0, kUnits::deltam32-kUnits::deltam21;
    return theReturn;
  }

  Eigen::Matrix3cd kIO(){// k matrix for mixing; normal ordering
    Eigen::Matrix3cd theReturn;
    theReturn << 0, 0, 0,
      0, kUnits::deltam21IO, 0,
      0, 0, kUnits::deltam32IO-kUnits::deltam21IO;
    return theReturn;
  }

  int factorial(int n){

    return (n==0) || (n==1) ? 1 : n* factorial(n-1);
  }
  
  Eigen::Matrix3cd sunSectionHamiltonian ( double minRad, double maxRad , double E){
    
    Eigen::Matrix3cd oneThing;// matrix for matter term
    oneThing << 1, 0, 0,
      0, 0, 0,
      0, 0, 0;
    //std::cout << "oneThing =\n" << oneThing << "\n";

    Eigen::Matrix3cd H_0 = U().transpose()*k()*U().conjugate(); // calculate U^* k U
    H_0 = H_0/(2*E); // finish prepping the vaccuum hamiltonian
  
    // fit from https://www.researchgate.net/publication/2225553_Solar_Models_Current_Epoch_and_Time_Dependences_Neutrinos_and_Helioseismological_Properties 
    // figure 8
    double expoConst = 245*kUnits::mol/kUnits::cm/kUnits::cm/kUnits::cm;
    double expoDecay = 10.54/kUnits::solarRadii;

    Eigen::Matrix3cd H_m = sqrt(2)*kUnits::fermiCouplingConstant*oneThing*(expoConst/expoDecay)*(exp(-expoDecay*minRad) - exp(-expoDecay*maxRad) )/(maxRad-minRad);
  
    Eigen::Matrix3cd H = H_0+H_m;// hamiltonian with matter term

    return H;
    
  }
  
  Eigen::Matrix3cd sunSectionHamiltonianIO ( double minRad, double maxRad , double E){
    
    Eigen::Matrix3cd oneThing;// matrix for matter term
    oneThing << 1, 0, 0,
      0, 0, 0,
      0, 0, 0;
    //std::cout << "oneThing =\n" << oneThing << "\n";

    Eigen::Matrix3cd H_0 = UIO().transpose()*kIO()*UIO().conjugate(); // calculate U^* k U
    H_0 = H_0/(2*E); // finish prepping the vaccuum hamiltonian
  
    // fit from https://www.researchgate.net/publication/2225553_Solar_Models_Current_Epoch_and_Time_Dependences_Neutrinos_and_Helioseismological_Properties 
    // figure 8
    double expoConst = 245*kUnits::mol/kUnits::cm/kUnits::cm/kUnits::cm;
    double expoDecay = 10.54/kUnits::solarRadii;

    Eigen::Matrix3cd H_m = sqrt(2)*kUnits::fermiCouplingConstant*oneThing*(expoConst/expoDecay)*(exp(-expoDecay*minRad) - exp(-expoDecay*maxRad) )/(maxRad-minRad);
  
    Eigen::Matrix3cd H = H_0+H_m;// hamiltonian with matter term

    return H;
    
  }

  Eigen::Matrix3cd myMatrixExp(Eigen::Matrix3cd theMatrix){
    double theMax = theMatrix.cwiseAbs().maxCoeff();
    long long int thePower = ceil(theMax);
    std::cout << "The power is " << thePower << ".\n";
    bool continuing = true;
    int iteration = 1;

    Eigen::Matrix3cd oldMat = Eigen::Matrix3cd::Identity();
    Eigen::Matrix3cd theMat = Eigen::Matrix3cd::Identity();
    do{
      
      Eigen::Matrix3cd tempMat = Eigen::Matrix3cd::Identity();
      for(int i = 0; i < iteration; i++){
	tempMat = tempMat*theMatrix/thePower;
      }
      theMat += tempMat/factorial(iteration);
      
      continuing = 2*(theMat.cwiseAbs().maxCoeff() - oldMat.cwiseAbs().maxCoeff())/(theMat.cwiseAbs().maxCoeff() + oldMat.cwiseAbs().maxCoeff()) > 0.001;
      
      oldMat = theMat;
      iteration++;
    }
    while(continuing);
    
    
    Eigen::Matrix3cd tempMat = Eigen::Matrix3cd::Identity();
    for(long long int i = 0; i < thePower; i++){
      tempMat = tempMat * theMat;
    }
    
    Eigen::Matrix3cd theReturn = tempMat;
    return theReturn;  
  }
  
  double* solarNeutrinoOscillationProbability(double E){
  
    
  
    double t = kUnits::solarRadii/kUnits::c; // time to escape the sun

    double n_e = 100*kUnits::mol/kUnits::cm/kUnits::cm/kUnits::cm;

    // initial H matrix; includes 1st term of e^x expansion
    Eigen::Matrix3cd timeEvolver = Eigen::Matrix3cd::Identity();
    //std::cout << "timeEvolver =\n" << timeEvolver << "\n before I start changing it.\n";


    Eigen::Matrix3cd oneThing;// matrix for matter term
    oneThing << 1, 0, 0,
      0, 0, 0,
      0, 0, 0;
    //std::cout << "oneThing =\n" << oneThing << "\n";

    Eigen::Matrix3cd H_0 = U().transpose()*k()*U().conjugate(); // calculate U^* k U
    H_0 = H_0/(2*E); // finish prepping the vaccuum hamiltonian
    //std::cout << "H_0 = U^* k U = \n" << H_0 << std::endl;
  

    Eigen::Matrix3cd H_m = sqrt(2)*kUnits::fermiCouplingConstant*n_e*oneThing;
  
    Eigen::Matrix3cd H = H_0+H_m;// hamiltonian with matter term
  

    timeEvolver = (-H*kUnits::i*t/kUnits::hbar).exp();
  
    Eigen::Vector3cd eKet;
    eKet << 1,0,0;
    Eigen::Vector3cd muKet;
    muKet << 0,1,0;
    Eigen::Vector3cd tauKet;
    tauKet << 0,0,1;
  
    std::complex<double> P_e = eKet.transpose() * timeEvolver * eKet;
    P_e = norm(P_e);
  
    std::complex<double> P_mu = muKet.transpose() * timeEvolver * eKet;
    P_mu = norm(P_mu);
  
    std::complex<double> P_tau = tauKet.transpose() * timeEvolver * eKet;
    P_tau = norm(P_tau);

    double *theReturn = new double[3];
    
    theReturn[0] = P_e.real();
    theReturn[1] = P_mu.real();
    theReturn[2] = P_tau.real();
    return theReturn;
  }
  
  double* solarNeutrinoSurvival(double E, double radius){
    bool pointFusion = false;// takes prescedence
    bool shellFusion = true;
    bool expoDensity = false;
    bool normalOrdering = true;
    int nSteps = 4;

    User::neutrinoData *theData = new User::neutrinoData;
    
    if(expoDensity){
      theData -> rStep = kUnits::solarRadii/nSteps;
      theData -> rNext = kUnits::solarRadii/nSteps;
      theData -> justMoved = 0;
    }
    
    if(pointFusion){
      theData -> neutrinoStartR = 0;
      theData -> neutrinoStartPhi = 0;
      theData -> neutrinoStartTheta = 0;
      
    }
    else if(shellFusion){
      double radius = 0.1*kUnits::solarRadii;// approximate radius for pp peak fusion
      User::shellFusionModel(radius, theData);
    }
    else{
      // put something here later for a relatively modular fusion model
    }
    
    theData -> X = theData -> neutrinoStartR*cos(theData -> neutrinoStartTheta)
      *sin(theData -> neutrinoStartPhi);
    //std::cout << "the starting X is " << (theData -> X)/kUnits::solarRadii << " solar radii.\n";
    //std::cin.ignore();
    theData -> Y = theData -> neutrinoStartR*sin(theData -> neutrinoStartTheta)
      *sin(theData -> neutrinoStartPhi);
    theData -> Z = theData -> neutrinoStartR*cos(theData -> neutrinoStartPhi);
    
    theData -> neutrinoEnergy = E;
    
    double n_e = 100*kUnits::mol/kUnits::cm/kUnits::cm/kUnits::cm;

    // initial H matrix; includes 1st term of e^x expansion
    Eigen::Matrix3cd timeEvolver = Eigen::Matrix3cd::Identity();


    Eigen::Matrix3cd oneThing;// matrix for matter term
    oneThing << 1, 0, 0,
      0, 0, 0,
      0, 0, 0;

    if(expoDensity){
      for(int i = 0; i < nSteps; i++){
	Eigen::Matrix3cd thisHamil;
	if(normalOrdering){
	  thisHamil = User::hamiltonianStep(theData);
	}
	else{
	  thisHamil = User::hamiltonianStepIO(theData);
	}
	
	timeEvolver = (-kUnits::i * thisHamil*(theData -> justMoved)/kUnits::c/kUnits::hbar).exp()
	  *timeEvolver; 
      }
      
      Eigen::Matrix3cd H_0;
      if(normalOrdering){
	H_0 = U().conjugate()*k()*U().transpose(); // calculate U^* k U^T
      }
      else{
	H_0 = UIO().conjugate()*kIO()*UIO().transpose(); // calculate U^* k U^T
      }

      H_0 = H_0/(2*E); // finish prepping the vaccuum hamiltonian
      
      timeEvolver = (-kUnits::i * H_0*(radius-1*kUnits::solarRadii)/kUnits::c/kUnits::hbar).exp() 
	* timeEvolver;
    }
    else{
      Eigen::Matrix3cd H_0;
      if(normalOrdering){
	H_0 = U().conjugate()*k()*U().transpose(); // calculate U^* k U
      }
      else{
	H_0 = UIO().conjugate()*kIO()*UIO().transpose(); // calculate U^* k U
      }
      H_0 = H_0/(2*E); // finish prepping the vaccuum hamiltonian
  
      Eigen::Matrix3cd H_m = sqrt(2)*kUnits::fermiCouplingConstant*n_e*oneThing;
      //std::cout << "The matter term is \n" << H_m << std::endl;
  
      timeEvolver = (-kUnits::i * (H_0 + H_m)*(kUnits::solarRadii + (theData -> Y) )/kUnits::c/kUnits::hbar).exp();

      timeEvolver = (-kUnits::i * H_0*(radius-1*kUnits::solarRadii)/kUnits::c/kUnits::hbar).exp() * timeEvolver;
    }
    delete theData;
  
    Eigen::Vector3cd eKet;
    eKet << 1,0,0;
    Eigen::Vector3cd muKet;
    muKet << 0,1,0;
    Eigen::Vector3cd tauKet;
    tauKet << 0,0,1;
  
    Eigen::Vector3cd eMassKet = U().conjugate().transpose()*eKet;
    Eigen::Vector3cd muMassKet = U().conjugate().transpose()*muKet;
    Eigen::Vector3cd tauMassKet = U().conjugate().transpose()*tauKet;

  
    std::complex<double> P_e = eKet.transpose() * timeEvolver * eKet;
    P_e = norm(P_e);
  
    std::complex<double> P_mu = muKet.transpose() * timeEvolver * eKet;
    P_mu = norm(P_mu);
  
    std::complex<double> P_tau = tauKet.transpose() * timeEvolver * eKet;
    P_tau = norm(P_tau);
    
    //std::cout << "P-e = \n" << P_e;
    //std::cin.ignore();

    double *theReturn = new double[3];
    
    //std::cout << "P(e)  = "<< P_e.real() << std::endl;
    //std::cout << "P(mu) = "<< P_e.real() << std::endl;
    //std::cout << "P(tau)= "<< P_e.real() << std::endl;
    
    theReturn[0] = P_e.real();
    theReturn[1] = P_mu.real();
    theReturn[2] = P_tau.real();
    return theReturn;
    
  }
}


/*
Eigen::Vector3cd* sunSurfacePacketMonteCarlo(double E, double Distance){
 
  Eigen::Matrix3cd H_0 = U().conjugate()*k()*U().transpose(); // calculate U^* k U
  H_0 = H_0/(2*E); // finish prepping the vaccuum hamiltonian
  
  Eigen::Matrix3cd H_m = sqrt(2)*kUnits::fermiCouplingConstant*n_e*oneThing;
  //std::cout << "The matter term is \n" << H_m << std::endl;
  
  Eigen::Vector3cd eKet;
  eKet << 1,0,0;
  Eigen::Vector3cd muKet;
  muKet << 0,1,0;
  Eigen::Vector3cd tauKet;
  tauKet << 0,0,1;

  Eigen::Matrix3cd matrixPart = (-kUnits::i * (H_0 + H_m)*(Distance )/kUnits::c/kUnits::hbar).exp()*eKet;
  
  double sigma_pP = kUnits::deltam32*kUnits::c*kUnits::c*kUnits::c/(2*E); // setting production momentum uncertainty to mass uncertainty scale from E^2=m^2+p^2
  
  double gaussConst = Pow(2*M_PI*sigma_pP^2,-1/4);

  // prepare integration
  Eigen::Matric3cd summation;// monte carlo integral: Sum at random points, divide by number of points, and multiply by the area of integration (min-max)
  summation << 0,0,0,0,0,0,0,0,0;
  std::default_random_engine generator;
  std::uniform_real_distribution<double> integralRange(E/kUnits::c-5*sigma_pP,E/kUnits::c+5*sigma_pP);// integrate to 5 st dev on gaussian





  long int i = 0;
  for(i = 0; i < 10; i++){
    summation += 
  }
  







  
  Eigen::Vector3cd* theReturn = new Eigen::Vector3cd;

  return theReturn;
  
}
*/

