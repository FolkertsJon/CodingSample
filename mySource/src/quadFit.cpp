

// This subprogram interpolates parabolas from the file we hand it
// We choose parabolas to account for the changing sign when our 
// orbit hits an edge. This function takes off two days.
// It also can return a vector containing points on an elliptical
// orbit following kepler's equation using newton's method

#include "quadFit.hh"

namespace User{
  
  // get number of lines  
  long long int number_of_lines (std::string theFile) {
    // char* intermediate = theFile;
    std::ifstream intermediate (theFile);
    //data.open (intermediate);
    std::cout << "I have opened the file \"" << theFile << "\"\n\n";
    std::string line;
    long long int nLines=0;
    while (std::getline(intermediate, line))
      ++nLines;
    std::cout << "Number of lines in text file: " << nLines << std::endl;
    return nLines;
  }
  
  double* johnTheInterpolator(std::string theFile){
    std::ifstream data;
    // get number of lines  
    long long int nLines = number_of_lines(theFile);

    double datum;
    data.open (theFile);
    std::cout << "I have opened the file \"" << theFile << " for the second time\"\n\n";
    long long int count=0;
    double y[nLines];// intermediate values  y(0 days), y(1 day), ...
    while(data >> datum){
      std::cout << "Line " << count << " gives datum " << datum << std::endl;
      y[count]=datum;
      count++;
    }
    std::cout << "I've collected the data set.\n";
    for(long long int i = 0; i<nLines; i++){
      //std::cout << intermediate[i] << endl;
    }

    // calculate a, b, c for at^2+bt+c in order to fit the radius as
    // a function of time for day 
    double *theValues = new double[3*(nLines-2)];// return the values as one long array
    for(long long int i=1; i<(nLines-1);i++){
      // values found from calculating (a,b,c)
      // =inverse( (1,x(i-1),x^2(i-1); 1,x(i),x^2(i)); 1,x(i+1),x^2(i+1) )
      // (inner product) (c,b,a)^T
    
      double a = y[i-1]/2-y[i]+y[i+1]/2;
      std::cout << "a[" << i << "] is " << a << std::endl;
      // ( 2*y[i]-y[i-1]-y[i+1] ) // numerator
      // / ( -2 ); // denominator
      double b = ( (-2*i-1)*y[i-1] + 4*i*y[i] + (-2*i+1)*y[i+1] )/2;
      /*( 2*y[i-1] - y[i] - y[i+1] ) / ( -3*i ) 
	- a * ( -6*i-1 ) / ( -3 );*/
      double c = ( (i*i+i)*y[i-1] + (2-2*i*i)*y[i] + (i*i-i)*y[i+1] )/2;
      //( y[i-1] + y[i] + y[i+1] + a*( 3*i+2 ) + b*( 3*i ) ) / 3;
      std::cout << "a[" << i << "] is " << a << ". b[" << i << "] is " << b << ". c[" << i << "] is " << c  << std::endl;
    
      // data for returning
      theValues [3*(i-1)]=c; // x^0 constant
      theValues [3*(i-1)+1]=b; // x^1 constant
      theValues [3*(i-1)+2]=a; // x^2 constant
    }

    return theValues;
  }


  std::vector<std::vector<double>> ellipticalRadiusTimeTest(double closest, double furthest){
    std::vector<std::vector<double>> radiusTime = {{},{}};
    std::cout << "I am " << furthest/kUnits::AU << " AU at my furthest.\n";
    std::cout << "I am " << closest/kUnits::solarRadii << " Rsol at my closest.\n";
 
    
    
    std::cout << "Furthest - closest = " << std::endl;
    std::cout << "Furthest + closest = " << std::endl;
    double e = (furthest-closest)/(furthest+closest);
    std::cout << "The ecentricity is " << e << ".\n";
    double a = (closest+furthest)/2;// semi-major axis
    std::cout << "The semi-major axis is is " << a << " solar radii.\n";
    double mu = -1;
    std::cout << "The mu parameter is " << mu << ".\n";
    
    
    // initial conditions
    double Radius = furthest;
    double Theta = -M_PI; //0
    double ThetaMax = M_PI;
    double dRadius = 0;
    // From Virial Theorem: 2T = -V
    double dTheta = 1e-6;//mu * (1+e*cos(Theta)) * (1+e*cos(Theta)) / (e * h * h * sin(Theta) * dRadius);

    // specific angular momentum; from def of angular momentum
    double h = 100;//dTheta * furthest * furthest;
    std::cout << "The h is " << h << ".\n";

    double dt = 1; // 
    double t = 0;

    double dRadiusMax = 0;

    while(Theta < ThetaMax){

      std::cout << Radius << ", " << dRadius << "/step, " << Theta << ", " << dTheta << "/step\n";
      // Save current values
      (radiusTime[0]).push_back(Radius);
      (radiusTime[1]).push_back(t);

      // update
      Radius += dRadius * dt;
      Theta += dTheta * dt;
      t += dt;
      std::cout << e * h * h * sin(Theta) * dTheta / ( (1 + e*cos(Theta)) * (1 + e*cos(Theta))) <<"\n";
      dRadius = e * h * h * sin(Theta) * dTheta / (mu * (1+ e*cos(Theta)) * (1 + e*cos(Theta)));
      dTheta = mu * (1+e*cos(Theta)) * (1+e*cos(Theta)) * dRadius / (e * h * h * sin(Theta));
      // save largest dRadius for checking if too large/small time steps
      if(dRadius > abs(dRadiusMax)) {
	dRadiusMax = dRadius;
      }
    }

    std::cout << "Total time  = " << t << "\n";
    
    
    return radiusTime;
  }



  std::vector<std::vector<double>> ellipticalRadiusTime(double closest, double furthest, double timeStep) {
    double a = (closest + furthest) / 2;
    double e = (furthest - closest) / (furthest + closest);

    double period = 2 * M_PI * sqrt(a * a * a / (kUnits::G * kUnits::mSun));
    //double timeStep = 100 * kUnits::sec;
    double endTime = period;
    double initialM = M_PI;

    std::vector<std::vector<double>> rtPairs;

    for (double t = 0; t <= endTime; t += timeStep) {
      double M = initialM + 2 * M_PI * t / period;
      double E = User::solveKepler(M, e);
      double r = User::computePositionForKepler(a, e, E);
      rtPairs.push_back({r, t});  // Store time and radius
    }

    return rtPairs;
  }

  


  double ellipticalRadius(double closest, double furthest, double time) {
    double a = (closest + furthest) / 2;
    double e = (furthest - closest) / (furthest + closest);

    double period = 2 * M_PI * sqrt(a * a * a / (kUnits::G * kUnits::mSun));

    double M = M_PI + 2 * M_PI * time / period; // add pi to start at furthest approach
    double E = User::solveKepler(M, e);
    double r = User::computePositionForKepler(a, e, E);

    return r;
  }

  /*
    std::vector<std::vector<double>> ellipticalRadiusTime(double closest, double furthest){
    std::vector<std::vector<double>> radiusTime = {{},{}};
    std::cout << "I am " << furthest/kUnits::AU << " AU at my furthest.\n";
    std::cout << "I am " << closest/kUnits::solarRadii << " Rsol at my closest.\n";

    
    
    std::cout << "Furthest - closest = " << ((furthest-closest)/kUnits::solarRadii) << std::endl;
    std::cout << "Furthest + closest = " << ((furthest+closest)/kUnits::solarRadii) << std::endl;
    double e = (furthest-closest)/(furthest+closest);
    std::cout << "The ecentricity is " << e << ".\n";
    double a = (closest+furthest)/2;// semi-major axis
    std::cout << "The semi-major axis is is " << a/kUnits::solarRadii << " solar radii.\n";
    double timeMax = 2*M_PI*sqrt( pow(a,3) / ( kUnits::G*kUnits::mSun ) );
    std::cout << "The period is " << timeMax/kUnits::day << " days.\n";
    double mu = -1*kUnits::G*kUnits::mSun;
    std::cout << "The mu parameter, -M_sun*G, is " << mu << " in SI units.\n";
    
    
    // initial conditions
    double Radius = closest;//furthest;
    double Theta = -M_PI; //0
    double dRadius = 0;
    // From Virial Theorem: 2T = -V
    double dTheta = sqrt(kUnits::G * kUnits::mSun/(furthest*furthest*furthest));

    // specific angular momentum; from def of angular momentum
    double h = abs(dTheta) * furthest * furthest;
    std::cout << "The h is " << h << ".\n";

    double dt = 10*kUnits::sec; // start at 10 second interaval; might change to variable time later
    double t = 0;

    size_t numSteps = round(timeMax/dt);

    double dRadiusMax = 0;

    for(size_t i = 0; i <= numSteps; i++){

    //std::cout << Radius/kUnits::solarRadii << ", " << dRadius*kUnits::sec/kUnits::solarRadii << ", " << Theta << ", " << dTheta*kUnits::sec << "\n";
    // Save current values
    (radiusTime[0]).push_back(Radius);
    (radiusTime[1]).push_back(t);

    // update
    Radius += dRadius * dt;
    Theta += dTheta * dt;
    t += dt;
    //std::cout << e * h * h * sin(Theta) * dTheta / ( (1 + e*cos(Theta)) * (1 + e*cos(Theta))) <<"\n";
    dRadius = e * h * h * sin(Theta) * dTheta / (mu * (1+ e*cos(Theta)) * (1 + e*cos(Theta)));
    dTheta = mu * (1+e*cos(Theta)) * (1+e*cos(Theta)) * dRadius / (e * h * h * sin(Theta));
    // save largest dRadius for checking if too large/small time steps
    if(dRadius*kUnits::sec/kUnits::solarRadii > abs(dRadiusMax)) {
    dRadiusMax = dRadius*kUnits::sec/kUnits::solarRadii;
    }
    }

    std::cout << "The largest dR/dt was : " << dRadiusMax << " Solar radii/sec\n";
    
    
    return radiusTime;
    }
  */

  double* elliptical (double closest, double furthest){
    std::cout << "I am " << furthest/kUnits::AU << " AU at my furthest.\n";
    std::cout << "I am " << closest/kUnits::solarRadii << " Rsol at my closest.\n";
  
    std::cout << "Furthest - closext = " << ((furthest-closest)/kUnits::solarRadii) << std::endl;
    std::cout << "Furthest + closext = " << ((furthest+closest)/kUnits::solarRadii) << std::endl;
    double e = (furthest-closest)/(furthest+closest);
    std::cout << "The eecentricity is " << e << ".\n";
    double p = furthest*(1-e);
    double a = (closest+furthest)/2;// semi-major axis
    double H = furthest*sqrt(kUnits::G*kUnits::mSun* ( 2/furthest - 1/a ) ) ;//  specific angular momentum
    const double timeStep = M_PI/10*kUnits::sec;
    std::cout << "The specific angular momentum is " << H << std::endl;

    double angle = -M_PI;// angle at the start; at furthest approach
    double radius =  p/(1 + e*cos(angle) );// r(theta)
    double timeMax = 2*M_PI*sqrt( pow(a,3) / ( kUnits::G*kUnits::mSun ) );
    std::cout << "\n\nI'm finding the orbital period by:  2*M_PI*sqrt( pow(a,3) / ( kUnits::G*kUnits::mSun ) ), where: \na = " << a << "\nG = " << kUnits::G << "\nmSun = " << kUnits::mSun << '\n';
    std::cout << "The orbital period is " << timeMax/kUnits::day << " days.\n";
    double nSteps = timeMax/timeStep;
    int intSteps = ceil(nSteps);
    long long int cells = ceil(timeMax/kUnits::min);
    std::cout << "The number of cells is " << cells << "\n";
    //double theValues[cells];
    double theValues[cells];
    int sanityCount = 0;


    std::cout << "I'm about to find the position/radii/etc for the ellipse uising " << intSteps << " steps. I am going to put them into " << cells << " cells.\n\n";
    for (long long int i = 0; i < intSteps; i++){
      int cellNum = ceil(i*timeStep/kUnits::min);    
      //    std::cout << "day = " << i*timeStep/kUnits::days << endl;
      if (i%10000==0){
	//std::cout << "r = " << radius/kUnits::solarRadii << std::endl;
	theValues[cellNum]=radius/kUnits::solarRadii;// in units of solar radii
	//std::cout << "I have written " << theValues[cellNum] << " to the array.\n";
      }
      double dAngle = H/(radius)/(radius);// dtheta/dt
      angle = angle + dAngle*timeStep;// update angle for dt
      radius = p/(1 + e*cos(angle) ); // update r for new angle



    
      //std::cout << "r = " << radius/kUnits::solarRadii << endl;
      //std::cout << "I'm filling cell number " << cellNum << " with " << radius/kUnits::solarRadii << endl;
      theValues[cellNum]=radius/kUnits::solarRadii;// in units of solar radii
    
      // Data updating
    
    }
    std::cout << "There are " << sizeof(theValues)/sizeof(theValues[0]) << " or " << cells << " values in the radius array.\n";

    // Check the radii
    for (int i=0;i<cells;i++){
      if(i%50000==0){
	std::cout << "The array value for the " << i << "th radius is r = " << theValues[i] << std::endl;
      }
    }


    std::cout << "\n\nI'm about to find the fit parameters for the ellipse, which should hold " << (3*(cells-1))  << " values.\n\n";
    double* theReturns = new double[(3*(cells-1))];// return the values as one long array
    // this type of declaration shuld fix a memory issue. We keep theReturns around until deleted


    for(long long int i=1; i<(cells-1);i++){
      // values found from calculating (a,b,c)
      // =inverse( (1,x(i-1),x^2(i-1); 1,x(i),x^2(i)); 1,x(i+1),x^2(i+1) )
      // (inner product) (c,b,a)^T
      double a = theValues[i-1]/2-theValues[i]+theValues[i+1]/2;
      // ( 2*theValues[i]-theValues[i-1]-theValues[i+1] ) // numerator
      // / ( -2 ); // denominator
      double b = ( (-2*i-1)*theValues[i-1] + 4*i*theValues[i] + (-2*i+1)*theValues[i+1] )/2;
      /*( 2*theValues[i-1] - theValues[i] - theValues[i+1] ) / ( -3*i ) 
	- a * ( -6*i-1 ) / ( -3 );*/
      double c = ( (i*i+i)*theValues[i-1] + (2-2*i*i)*theValues[i] + (i*i-i)*theValues[i+1] )/2;
    
      if( i%100==0){
	//std::cout << "c = " << c <<".\n";
      }

      if (c>1000000){
	std::cout << "c is more than a million... somehow. The iteration number is "
		  << i << ", and the cell number is " << (3*(i-1)) << ".\nFor good measure, b = " 
		  << b <<", and a = " << a << std::endl;
	break;
      }  
      if (c<-1000000){
	std::cout << "c is less than negative 1 million... somehow. The cell number is "	   << i << ".\n";
	break;
      }  
      //( theValues[i-1] + theValues[i] + theValues[i+1] + a*( 3*i+2 ) + b*( 3*i ) ) / 3;
      //std::cout << "a[" << i-1 << "] is " << a << ". b[" << i-1 << "] is " << b << ". c[" << i-1 << "] is " << c  << endl;
      //std::cout << "Hence the radius at day " << i << " is " << (a*i*i + b*i + c    ) << ".\n\n";

      // data for returning
      theReturns [3*(i-1)]=c; // x^0 constant
      theReturns [3*(i-1)+1]=b; // x^1 constant
      theReturns [3*(i-1)+2]=a; // x^2 constant
    }
    std::cout << "c[0], b[0], and a[0], are " << theReturns[0] << ", " << theReturns[1] << ", " << theReturns[2] << " respectively.\n\n";

    std::cout << "elliptical is done iterating.\n";
  
    return theReturns;
  }




  double positionDays(double time, double* quadConst, double maxTime){
    //std::cout << "Position function has been called.\n";
    double currentRadius=0;
    int day = floor(time/kUnits::day);
    /* std::cout is slow
       std::cout << "The time is " << (time/kUnits::hr) << "hours in, and we are " 
       << day << " days in.\n";
    */
    if(time < (1*kUnits::day)){// first day from later data
      currentRadius = *(quadConst)+time*(*( quadConst+1 ) )/kUnits::day 
	+ time * time * (*( quadConst+2 ) ) /kUnits::day /kUnits::day;
    }
    else if (time > (maxTime-kUnits::day)){ // last day from earlier data
      currentRadius = *(quadConst + 3*( day-2 ) )
	+ *( quadConst+1  + 3*( day-2 ) ) * time /kUnits::day 
	+ *( quadConst+2  + 3*( day-2 ) ) * time * time /kUnits::day /kUnits::day;    
    }
    else{
      currentRadius = *(quadConst + 3*( day-1 ) )
	+ *( quadConst+1  + 3*( day-1 ) ) * time /kUnits::day 
	+ *( quadConst+2  + 3*( day-1 ) ) * time * time /kUnits::day /kUnits::day;
    }
    //std::cout << "The current radius is " << currentRadius << endl;
    return currentRadius*kUnits::solarRadii; // return the radius in solar radii
  }

  // hours version
  double positionHours(double time, double* quadConst, double maxTime){
    //std::cout << "Position function has been called.\n";
    double currentRadius=0;
    int day = floor(time/kUnits::day);
    //std::cout << "The time is " << (time/kUnits::hr) << "hours in, and we are " 
    //     << day << " days in.\n";
    if(time < (1*kUnits::hr)){
      currentRadius = *(quadConst)+time*(*( quadConst+1 ) )/kUnits::hr 
	+ time * time * (*( quadConst+2 ) ) /kUnits::hr /kUnits::hr;
    }
    else if (time > (maxTime - kUnits::hr)){ // last day from earlier data
      currentRadius = *(quadConst + 3*( day-2 ) )
	+ *( quadConst+1  + 3*( day-2 ) ) * time /kUnits::hr 
	+ *( quadConst+2  + 3*( day-2 ) ) * time * time /kUnits::hr /kUnits::hr;    
      /*      
	      std::cout << "This is the last time set. a = " 
	      << *( quadConst+2  + 3*( day-2 ) ) << ". b = " 
	      << *( quadConst+1  + 3*( day-2 ) ) << ". c = " 
	      << *( quadConst+2  + 3*( day-2 ) ) << ". t = " << time 
	      << ". Hence the radius is " 
	      << (*(quadConst + 3*( day-2 ) )
	      + *( quadConst+1  + 3*( day-2 ) ) * time /kUnits::hr 
	      + *( quadConst+2  + 3*( day-2 ) ) * time * time /kUnits::hr /kUnits::hr)
	      << ".\n";
      */
    }
    else{
      int day = floor(time/kUnits::hr);
      /* std::couts are slow
	 std::cout << "a = " << *( quadConst+2  + 3*( day-1 ) ) << ". b = "
	 << *( quadConst+1  + 3*( day-1 ) ) << ". c = "
	 << *(quadConst + 3*( day-1 ) ) << ".\n";
      */

      currentRadius = *(quadConst + 3*( day-1 ) )
	+ *( quadConst+1  + 3*( day-1 ) ) * time /kUnits::hr 
	+ *( quadConst+2  + 3*( day-1 ) ) * time * time /kUnits::hr /kUnits::hr;
    
    }
    //std::cout << "The current radius is " << currentRadius << endl;
    return currentRadius*kUnits::solarRadii; // return the radius in solar radii
  }

  // minutes version
  double positionMinutes(double time, double* quadConst, double maxTime){
    //std::cout << "Position function has been called.\n";
    double currentRadius=0;
    int day = floor(time/kUnits::day);
    //std::cout << "The time is " << (time/kUnits::hr) << "hours in, and we are " 
    //     << day << " days in.\n";
    if(time < (1*kUnits::min)){
      currentRadius = *(quadConst)+time*(*( quadConst+1 ) )/kUnits::min 
	+ time * time * (*( quadConst+2 ) ) /kUnits::min /kUnits::min;
    }
    else if (time > (maxTime - kUnits::hr)){ // last day from earlier data
      currentRadius = *(quadConst + 3*( day-2 ) )
	+ *( quadConst+1  + 3*( day-2 ) ) * time /kUnits::minute 
	+ *( quadConst+2  + 3*( day-2 ) ) * time * time /kUnits::minute /kUnits::min;    
    
      /*
	std::cout << "This is the last time set. a = " << *( quadConst+2  + 3*( day-2 ) ) 
	<< ". b = " << *( quadConst+1  + 3*( day-2 ) ) << ". c = " 
	<< *( quadConst+2  + 3*( day-2 ) ) << ". t = " << time 
	<< ". Hence the radius is " << (*(quadConst + 3*( day-2 ) )
	+ *( quadConst+1  + 3*( day-2 ) ) * time /kUnits::minute 
	+ *( quadConst+2  + 3*( day-2 ) ) * time * time /kUnits::minute /kUnits::min)
	<< ".\n";
      */
    }
    else{
      int day = floor(time/kUnits::min);
      /* std::couts are slow
	 std::cout << "a = " << *( quadConst+2  + 3*( day-1 ) ) << ". b = "
	 << *( quadConst+1  + 3*( day-1 ) ) << ". c = "
	 << *(quadConst + 3*( day-1 ) ) << ".\n";
      */
      currentRadius = *(quadConst + 3*( day-1 ) )
	+ *( quadConst+1  + 3*( day-1 ) ) * time /kUnits::min 
	+ *( quadConst+2  + 3*( day-1 ) ) * time * time /kUnits::min /kUnits::min;
    
    }
    // std::cout << "The current radius is " << currentRadius << endl;
    return currentRadius*kUnits::solarRadii; // return the radius in solar radii
  }

  double derivativeMinutes(double time, double* quadConst, double maxTime){
    double derivative = 0;
    int day = floor(time/kUnits::day);
    if(time < (1*kUnits::min)){
      derivative = (*( quadConst+1 ) ) 
	+ 2 * time * (*( quadConst+2 ) )  /kUnits::min;
    }
    else if (time > (maxTime - kUnits::hr)){ // last day from earlier data
      derivative = 
	+ *( quadConst+1  + 3*( day-2 ) ) / kUnits::minute 
	+ *( quadConst+2  + 3*( day-2 ) ) * 2 * time /kUnits::min;    
    }
    else{
      int day = floor(time/kUnits::min);
      derivative = 
	+ *( quadConst+1  + 3*( day-1 ) ) / kUnits::min 
	+ *( quadConst+2  + 3*( day-1 ) ) * 2 * time /kUnits::min;
    
    }
    if(derivative < 0) derivative = -derivative;
    return derivative*kUnits::solarRadii; // return the deravative in solar radii
  }



  std::vector<std::pair<double, double>> radiusTimePairsFromFile(std::string filename){
    std::vector<std::pair<double, double>> pairs;
    std::ifstream file(filename.c_str());
    std::string line;
    double t, r;

    if (!file.is_open()) {
      std::cerr << "Failed to open file: " << filename << std::endl;
      return pairs; // Return empty vector if file cannot be opened
    }

    while (std::getline(file, line)) {
      std::istringstream iss(line);
      char comma; // to read the comma separator

      if (!(iss >> t >> comma >> r)) {
	std::cerr << "Failed to parse line: " << line << std::endl;
	continue; // Skip the malformed line
      }

      pairs.push_back(std::make_pair(r*kUnits::AU, t*kUnits::day));
    }

    file.close();
    return pairs;
  }



  double interpPosition(std::vector<std::pair<double, double>> pairs, double timeToInterpolate) {
    if (pairs.empty()) {
      std::cerr << "The vector of pairs is empty." << std::endl;
      return 0.0;
    }
    if (pairs.size() == 1) {
      return pairs[0].first; // Return the radius of the single element
    }

    for (size_t i = 0; i < pairs.size() - 1; i++) {
      if (timeToInterpolate >= pairs[i].second && timeToInterpolate <= pairs[i + 1].second) {
	double t1 = pairs[i].second, r1 = pairs[i].first;
	double t2 = pairs[i + 1].second, r2 = pairs[i + 1].first;
	return r1 + (r2 - r1) * (timeToInterpolate - t1) / (t2 - t1); // Linear interpolation formula
      }
    }

    // Extrapolation part
    // If the time to interpolate is after the last time point, extrapolate using the last two points
    double t1 = pairs[pairs.size() - 2].second, r1 = pairs[pairs.size() - 2].first;
    double t2 = pairs[pairs.size() - 1].second, r2 = pairs[pairs.size() - 1].first;
    if (timeToInterpolate > t2) {
      // Extrapolate based on the slope of the last two points
      return r2 + (r2 - r1) * (timeToInterpolate - t2) / (t2 - t1);
    }

    // If the specified time is before the first time point, handle as desired
    if (timeToInterpolate < pairs[0].second) {
      std::cerr << "Time to interpolate (" << timeToInterpolate << ") is before the range of the provided data." << std::endl;
      return 0.0; // Returning zero or any default value you think is appropriate
    }

    return 0.0; // Default return in case all else fails
  }





  
  
}// End User Namespace
