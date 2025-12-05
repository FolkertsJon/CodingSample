

// This subprogram interpolates parabolas from the file we hand it
// We choose parabolas to account for the changing sign when our 
// orbit hits an edge. This function takes off two days.

#include "kUnits.hh"
#include "quadraticFits.hh"


double *johnTheInterpolator(std::string theFile){
  ifstream data;
  // get number of lines  
  int nLines = number_of_lines(theFile);

  double datum;
  data.open (theFile);
  cout << "I have opened the file \"" << theFile << " for the second time\"\n\n";
  int count=0;
  double y[nLines];// intermediate values  y(0 days), y(1 day), ...
  while(data >> datum){
    cout << "Line " << count << " gives datum " << datum <<endl;
    y[count]=datum;
    count++;
  }
  cout << "I've collected the data set.\n";
  for(int i = 0; i<nLines; i++){
    //cout << intermediate[i] << endl;
  }

  // calculate a, b, c for at^2+bt+c in order to fit the radius as
  // a function of time for day 
  double *theValues = new  double[3*(nLines-2)];// return the values as one long array
  for(int i=1; i<(nLines-1);i++){
    // values found from calculating (a,b,c)
    // =inverse( (1,x(i-1),x^2(i-1); 1,x(i),x^2(i)); 1,x(i+1),x^2(i+1) )
    // (inner product) (c,b,a)^T
    
    double a = y[i-1]/2-y[i]+y[i+1]/2;
    cout << "a[" << i << "] is " << a << endl;
      // ( 2*y[i]-y[i-1]-y[i+1] ) // numerator
      // / ( -2 ); // denominator
    double b = ( (-2*i-1)*y[i-1] + 4*i*y[i] + (-2*i+1)*y[i+1] )/2;
      /*( 2*y[i-1] - y[i] - y[i+1] ) / ( -3*i ) 
	- a * ( -6*i-1 ) / ( -3 );*/
    double c = ( (i*i+i)*y[i-1] + (2-2*i*i)*y[i] + (i*i-i)*y[i+1] )/2;
      //( y[i-1] + y[i] + y[i+1] + a*( 3*i+2 ) + b*( 3*i ) ) / 3;
    cout << "a[" << i << "] is " << a << ". b[" << i << "] is " << b << ". c[" << i << "] is " << c  << endl;
    
    // data for returning
    theValues [3*(i-1)]=c; // x^0 constant
    theValues [3*(i-1)+1]=b; // x^1 constant
    theValues [3*(i-1)+2]=a; // x^2 constant
  }

  return theValues;
}



double* elliptical (double closest, double furthest){
  cout << "I am " << furthest/kUnits::AU << " AU at my furthest.\n";
  cout << "I am " << closest/kUnits::solarRadii << " Rsol at my closest.\n";
  
  cout << "Furthest - closext = " << ((furthest-closest)/kUnits::solarRadii) << endl;
  cout << "Furthest + closext = " << ((furthest+closest)/kUnits::solarRadii) << endl;
  double e = (furthest-closest)/(furthest+closest);
  cout << "The eecentricity is " << e << ".\n";
  double p = furthest*(1-e);
  double a = (closest+furthest)/2;// semi-major axis
  double H = furthest*sqrt(kUnits::G*kUnits::mSun* ( 2/furthest - 1/a ) ) ;//  specific angular momentum
  const double timeStep = M_PI/10*kUnits::sec;  cout << "The specific angular momentum is " << H << endl;

  double angle = -M_PI;// angle at the start; at furthest approach
  double radius =  p/(1 + e*cos(angle) );// r(theta)
  double timeMax = 2*M_PI*sqrt( pow(a,3) / ( kUnits::G*kUnits::mSun ) );
  cout << "The orbital period is " << timeMax/kUnits::day << " days.\n";
  double nSteps = timeMax/timeStep;
  int intSteps = ceil(nSteps);
  long long int cells = ceil(timeMax/kUnits::min);
  //double theValues[cells];
  double theValues[cells];
  int sanityCount = 0;


  cout << "I'm about to find the position/radii/etc for the ellipse uising " << intSteps << " steps. I am going to put them into " << cells << " cells.\n\n";
  for (long long int i = 0; i < intSteps; i++){
    int cellNum = ceil(i*timeStep/kUnits::min);    
    //    cout << "day = " << i*timeStep/kUnits::days << endl;
    if (i%10000==0){
    cout << "r = " << radius/kUnits::solarRadii << endl;
    theValues[cellNum]=radius/kUnits::solarRadii;// in units of solar radii
    //cout << "I have written " << theValues[cellNum] << " to the array.\n";
    }
    double dAngle = H/(radius)/(radius);// dtheta/dt
    angle = angle + dAngle*timeStep;// update angle for dt
    radius = p/(1 + e*cos(angle) ); // update r for new angle



    
    //cout << "r = " << radius/kUnits::solarRadii << endl;
    //cout << "I'm filling cell number " << cellNum << " with " << radius/kUnits::solarRadii << endl;
    theValues[cellNum]=radius/kUnits::solarRadii;// in units of solar radii
    
    // Data updating
    
  }
  cout << "There are " << sizeof(theValues)/sizeof(theValues[0]) << " or " << cells << " values in the radius array.\n";

  // Check the radii
  for (int i=0;i<cells;i++){
    if(i%50000==0){
      cout << "The array value for the " << i << "th radius is r = " << theValues[i] << endl;
    }
  }


  cout << "\n\nI'm about to find the fit parameters for the ellipse, which should hold " << (3*(cells-1))  << "values.\n\n";
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
    cout << "c = " << c <<".\n";
    }

    if (c>1000000){
      cout << "c is more than a million... somehow. The iteration number is "
	   << i << ", and the cell number is " << (3*(i-1)) << ".\nFor good measure, b = " 
	   << b <<", and a = " << a << endl;
	break;
    }  
    if (c<-1000000){
      cout << "c is less than negative 1 million... somehow. The cell number is "	   << i << ".\n";
	break;
    }  
    //( theValues[i-1] + theValues[i] + theValues[i+1] + a*( 3*i+2 ) + b*( 3*i ) ) / 3;
    //cout << "a[" << i-1 << "] is " << a << ". b[" << i-1 << "] is " << b << ". c[" << i-1 << "] is " << c  << endl;
    //cout << "Hence the radius at day " << i << " is " << (a*i*i + b*i + c    ) << ".\n\n";

    // data for returning
    theReturns [3*(i-1)]=c; // x^0 constant
    theReturns [3*(i-1)+1]=b; // x^1 constant
    theReturns [3*(i-1)+2]=a; // x^2 constant
  }
  cout << "c[0], b[0], and a[0], are " << theReturns[0] << ", " << theReturns[1] << ", " << theReturns[2] << " respectively.\n\n";

  cout << "elliptical is done iterating.\n";
  
  return theReturns;
}




double positionDays(double time, double* quadConst, double maxTime){
  //cout << "Position function has been called.\n";
  double currentRadius=0;
  int day = floor(time/kUnits::day);
  /* cout is slow
    cout << "The time is " << (time/kUnits::hr) << "hours in, and we are " 
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
  //cout << "The current radius is " << currentRadius << endl;
  return currentRadius*kUnits::solarRadii; // return the radius in solar radii
}

// hours version
double positionHours(double time, double* quadConst, double maxTime){
  //cout << "Position function has been called.\n";
  double currentRadius=0;
  int day = floor(time/kUnits::day);
  //cout << "The time is " << (time/kUnits::hr) << "hours in, and we are " 
  //     << day << " days in.\n";
  if(time < (1*kUnits::hr)){
    currentRadius = *(quadConst)+time*(*( quadConst+1 ) )/kUnits::hr 
      + time * time * (*( quadConst+2 ) ) /kUnits::hr /kUnits::hr;
  }
  else if (time > (maxTime - kUnits::hr)){ // last day from earlier data
    currentRadius = *(quadConst + 3*( day-2 ) )
      + *( quadConst+1  + 3*( day-2 ) ) * time /kUnits::hr 
      + *( quadConst+2  + 3*( day-2 ) ) * time * time /kUnits::hr /kUnits::hr;    
    cout << "This is the last time set. a = " 
	 << *( quadConst+2  + 3*( day-2 ) ) << ". b = " 
	 << *( quadConst+1  + 3*( day-2 ) ) << ". c = " 
	 << *( quadConst+2  + 3*( day-2 ) ) << ". t = " << time 
	 << ". Hence the radius is " 
	 << (*(quadConst + 3*( day-2 ) )
	     + *( quadConst+1  + 3*( day-2 ) ) * time /kUnits::hr 
	     + *( quadConst+2  + 3*( day-2 ) ) * time * time /kUnits::hr /kUnits::hr)
	 << ".\n";
  }
  else{
    int day = floor(time/kUnits::hr);
    /* couts are slow
    cout << "a = " << *( quadConst+2  + 3*( day-1 ) ) << ". b = "
	 << *( quadConst+1  + 3*( day-1 ) ) << ". c = "
	 << *(quadConst + 3*( day-1 ) ) << ".\n";
    */

    currentRadius = *(quadConst + 3*( day-1 ) )
      + *( quadConst+1  + 3*( day-1 ) ) * time /kUnits::hr 
      + *( quadConst+2  + 3*( day-1 ) ) * time * time /kUnits::hr /kUnits::hr;
    
  }
  //cout << "The current radius is " << currentRadius << endl;
  return currentRadius*kUnits::solarRadii; // return the radius in solar radii
}

// minutes version
double positionMinutes(double time, double* quadConst, double maxTime){
  //cout << "Position function has been called.\n";
  double currentRadius=0;
  int day = floor(time/kUnits::day);
  //cout << "The time is " << (time/kUnits::hr) << "hours in, and we are " 
  //     << day << " days in.\n";
  if(time < (1*kUnits::min)){
    currentRadius = *(quadConst)+time*(*( quadConst+1 ) )/kUnits::min 
      + time * time * (*( quadConst+2 ) ) /kUnits::min /kUnits::min;
  }
  else if (time > (maxTime - kUnits::hr)){ // last day from earlier data
    currentRadius = *(quadConst + 3*( day-2 ) )
      + *( quadConst+1  + 3*( day-2 ) ) * time /kUnits::minute 
      + *( quadConst+2  + 3*( day-2 ) ) * time * time /kUnits::minute /kUnits::min;    
    
    cout << "This is the last time set. a = " << *( quadConst+2  + 3*( day-2 ) ) 
	 << ". b = " << *( quadConst+1  + 3*( day-2 ) ) << ". c = " 
	 << *( quadConst+2  + 3*( day-2 ) ) << ". t = " << time 
	 << ". Hence the radius is " << (*(quadConst + 3*( day-2 ) )
					 + *( quadConst+1  + 3*( day-2 ) ) * time /kUnits::minute 
					 + *( quadConst+2  + 3*( day-2 ) ) * time * time /kUnits::minute /kUnits::min)
	 << ".\n";
  }
  else{
    int day = floor(time/kUnits::min);
    /* couts are slow
      cout << "a = " << *( quadConst+2  + 3*( day-1 ) ) << ". b = "
           << *( quadConst+1  + 3*( day-1 ) ) << ". c = "
	   << *(quadConst + 3*( day-1 ) ) << ".\n";
    */
    currentRadius = *(quadConst + 3*( day-1 ) )
      + *( quadConst+1  + 3*( day-1 ) ) * time /kUnits::min 
      + *( quadConst+2  + 3*( day-1 ) ) * time * time /kUnits::min /kUnits::min;
    
  }
  // cout << "The current radius is " << currentRadius << endl;
  return currentRadius*kUnits::solarRadii; // return the radius in solar radii
}
