// This file contains all the functions to take in data from the veto
// SiPMs used for the LS and veto testing during summer 2022. The
// data comes from csv files output by the digital oscilloscope, but
// the data are paresed before being passed to this 
//needed to call from main.cc
// in order to change the parameters (time, filemode, etc)
// This contains the functions that can check for arguments, parse
// the argument(s) that are input, and then execute the appropriate
// changes.
// Jonathan Folkerts: Created June 2022

#include <pulseShaper.hh>

namespace User{

  // This function will take in a vector with the pulse, and output 
  // the shaped pulse
  std::vector<double> *vetoPulseShaper( std::vector<double>* pulseIn, std::vector<double>* pulseTime){
    int vecSize = pulseIn->size();
    // Don't iterate through twice
    double maximum = 0; //*max_element( pulseIn->begin(), pulseIn->end() );
    double maxTime = -1*kUnits::sec;
    int maxIdx = vecSize;
    

    double timeStep = 0*kUnits::sec;
    if( double k = (*pulseTime)[1] ){
      timeStep = (*pulseTime)[1] - (*pulseTime)[0];
    }

    double eVal = 0; // Maximum/euler's constant
    double eTime = -1*kUnits::sec;
    int eIdx = vecSize;

    // vector to hold derivative of pulse shape
    std::vector<double> *derivative = new std::vector<double>(vecSize);
    std::vector<double> *avgDev = new std::vector<double>(vecSize);
    
    // Find and set the maximum value, index, and time
    for(int i = 0; i < vecSize; i++){
      // Find the derivative
      if( (i > 0) && ( i < ( vecSize - 1 ) ) ){
	// central difference
	(*derivative)[i] = ( (*pulseIn)[i+1] - (*pulseIn)[i-1] )/(2*timeStep);
      }
      else if( i = 0 ){
	// forward difference
	(*derivative)[i] = ( (*pulseIn)[i+1] - (*pulseIn)[i] )/(timeStep);
      }
      else if( i = (vecSize - 1) ){
	// backward difference
	(*derivative)[i] = ( (*pulseIn)[i] - (*pulseIn)[i-1] )/(timeStep);
      }
      
      
      // Find the max
      if( (*pulseIn)[i] > maximum ){
	maximum = (*pulseIn)[i];
	maxIdx = i;
      }
    }
    maxTime = (*pulseTime)[maxIdx]*kUnits::sec;


    // Find and set the 1/e value, time and index
    eVal = maximum*exp(-1);
    double dummy = 0;
    for(int i=0; i < vecSize; i++){
      // Find the averaged derivative
      if( (i > 0) && ( i < ( vecSize - 1 ) ) ){
	// central 
	(*avgDev)[i] = ( (*pulseIn)[i+1] + 2*(*pulseIn)[i] + (*pulseIn)[i-1] )/4;
      }
      else if( i = 0 ){
	// forward
	(*avgDev)[i] = ( (*pulseIn)[i+1] + 2*(*pulseIn)[i] )/3;
      }
      else if( i = (vecSize - 1) ){
	// backward 
	(*avgDev)[i] = ( 2*(*pulseIn)[i] - (*pulseIn)[i-1] )/3;
      }
      
      if( ( (*pulseIn)[i] > dummy) &&  ( (*pulseIn)[i] < eVal) && ( i < maxIdx ) ){
	dummy = (*pulseIn)[i];
	eIdx = i;
      }
    }
    eTime = (*pulseTime)[eIdx]*kUnits::sec;
    
    double leadTime = 20*kUnits::ns;
    double muPulseIntegral = 100*kUnits::mV*kUnits::ns;
    
    
    
    
    
    
    
    
    // Cleanup
    delete derivative;
    delete avgDev;
    // output
    std::vector<double> test = {1,2,3};
    std::vector<double> *output = new std::vector<double>(test);
  }







std::vector<double> *signalMaxTimeAndMax( std::vector<double>* pulseIn, std::vector<double>* pulseTime){
    int vecSize = pulseIn->size();
    //std::cout << "the vector size is " << vecSize << "\n";
    // Don't iterate through twice
    


    double maximum = 16e-4; //*max_element( pulseIn->begin(), pulseIn->end() );
    double maxTime = -1*kUnits::sec;
    int maxIdx = 0;
    

    double timeStep = 0*kUnits::sec;
    if( double k = (*pulseTime)[1] ){
      timeStep = (*pulseTime)[1] - (*pulseTime)[0];
    }

    
    // Find and set the maximum value, index, and time
    for(int i = 0; i < vecSize; i++){      
      // Find the max
	
      if( i % 1000 == 0 ){
	//std::cout << "(*pulseIn)[" <<i << "] = " << (*pulseIn)[i] << "\n";
      }

      if( (*pulseIn)[i] > maximum ){
	//std::cout << "The maximu is being set to " << (*pulseIn)[i] << "\n";
	maximum = (*pulseIn)[i];
	maxIdx = i;
      }
    }
    maxTime = (*pulseTime)[maxIdx]*kUnits::sec;


    
    
    
    
    //std::cout << "1/e -> max time = " << RCTime << ", and maximum value = " << maximum << "\n";
    
    // output
    std::vector<double> initializer = {maxTime,maximum};
    std::vector<double> *output = new std::vector<double>(initializer);
    return output;
  }





std::vector<double> *vetoRCTimeAndMax( std::vector<double>* pulseIn, std::vector<double>* pulseTime, std::vector<bool>* vetoArr){
    int vecSize = pulseIn->size();
    //std::cout << "the vector size is " << vecSize << "\n";
    // Don't iterate through twice
    


    double maximum = 16e-4; //*max_element( pulseIn->begin(), pulseIn->end() );
    double maxTime = -1*kUnits::sec;
    int maxIdx = 0;
    

    double timeStep = 0*kUnits::sec;
    if( double k = (*pulseTime)[1] ){
      timeStep = (*pulseTime)[1] - (*pulseTime)[0];
    }
    
    double vetoTime = 150*kUnits::ns; // length of veto pulse
    int vetoCountMax = floor(vetoTime/timeStep);
    int vetoCount = 0;

    double offsetTime = 20*kUnits::ns; // the offest we're going to move the veto earlier by
    int offsetCount = floor(offsetTime/timeStep);

    double eVal = 16e-4; // Maximum/euler's constant
    double eTime = -1*kUnits::sec;
    int eIdx = 0;
    
    // Find and set the maximum value, index, and time
    for(int i = offsetCount; i < (vecSize-offsetCount); i++){ // skip the first and last offsetCount indices
      // Find the max
	
      if( i % 1000 == 0 ){
	//std::cout << "(*pulseIn)[" <<i << "] = " << (*pulseIn)[i] << "\n";
      }
      
      bool veto = (*vetoArr)[i];
      if( veto && vetoCount == 0){// check to see if we're vetoing
	vetoCount = vetoCountMax;
      }
      if( vetoCount == 0 ){ // only if veto is false
	if( (*pulseIn)[i] > maximum ){
	  //std::cout << "The maximu is being set to " << (*pulseIn)[i] << "\n";
	  maximum = (*pulseIn)[i];
	  maxIdx = i;
	}
      }
      vetoCount--;
    }
    maxTime = (*pulseTime)[maxIdx]*kUnits::sec;


    // Find and set the 1/e value, time and index
    eVal = maximum*exp(-1);
    double dummy = 0;
    for(int i = maxIdx; i > 0 ; i--){// searching backwards from the maxIdx
      // Find the 1/e value
      if( ( (*pulseIn)[i] > dummy) &&  ( (*pulseIn)[i] < eVal) && ( i < maxIdx ) ){
	//std::cout << "The 1/e is being set to " << (*pulseIn)[i] << " Which corresponds to time " << (*pulseTime)[i] << "\n";
	dummy = (*pulseIn)[i];
	eIdx = i;
      }
    }
    eTime = (*pulseTime)[eIdx]*kUnits::sec;
    
    double RCTime = maxTime - eTime;
    
    
    
    
    //std::cout << "1/e -> max time = " << RCTime << ", and maximum value = " << maximum << "\n";
    
    // output
    std::vector<double> test = {RCTime,maximum};
    std::vector<double> *output = new std::vector<double>(test);
    return output;
  }





std::vector<double> *vetoRCTimeAndMax( std::vector<double>* pulseIn, std::vector<double>* pulseTime){
  int size = pulseIn -> size();
  std::vector<bool> temp(size,false);
  return vetoRCTimeAndMax( pulseIn, pulseTime, &temp);
  }





  
  std::vector<bool> *produceVeto( std::vector<double>* vetoPulse, double timeStep){
    int vecSize = vetoPulse -> size();
    bool vetoing = false;
    double threshold = 0.01*kUnits::V; // value subject to change
    int countdown = 0; 

    std::vector<bool> boolArray(vecSize, false);
    
    for( int i = 0; i < vecSize; i++){
      // check for and set vetoing bool
      if( (*vetoPulse)[i] > threshold){
	vetoing = true;
      }
      else{
	vetoing = false;
      }
      
      // reset the count if we're still in vetoland
      if(vetoing && (countdown == 0)){
	countdown = floor(200*kUnits::ns/timeStep);
      }
      
      // set array value true and count down
      if( countdown > 0){
	boolArray[i] = true;
	
	countdown--;
      }
      
    }
    
    std::vector<bool> *theReturn = new std::vector<bool>(boolArray);
    
  }





  double pulseArea(std::vector<double>* pulseIn, double stepSize){
    int vecSize = pulseIn->size();

    double maximum = 16e-4; //*max_element( pulseIn->begin(), pulseIn->end() );
    int maxIdx = 0;
    
    double eVal = 16e-4; // Maximum/euler's constant
    int eIdx1 = 0;
    int eIdx2 = 0;
    
    // Find and set the maximum value and index
    for(int i = 0; i < vecSize; i++){      
      // Find the max
	
      if( i % 1000 == 0 ){
	//std::cout << "(*pulseIn)[" <<i << "] = " << (*pulseIn)[i] << "\n";
      }

      if( (*pulseIn)[i] > maximum ){
	//std::cout << "The maximu is being set to " << (*pulseIn)[i] << "\n";
	maximum = (*pulseIn)[i];
	maxIdx = i;
      }
    }
    
    
    // Find and set the 1/e value and index
    eVal = maximum*exp(-1);
    double dummy = 0;
    for(int i = maxIdx; i > 0 ; i--){// searching backwards from the maxIdx
      // Find the 1/e value
      if( ( (*pulseIn)[i] > dummy) &&  ( (*pulseIn)[i] < eVal) && ( i < maxIdx ) ){
	dummy = (*pulseIn)[i];
	eIdx1 = i;
      }
    }
    dummy = 0;
    for(int i = maxIdx; i < vecSize ; i++){// searching forwards from the maxIdx
      // Find the 1/e value
      if( ( (*pulseIn)[i] > dummy) &&  ( (*pulseIn)[i] < eVal) && ( i > maxIdx ) ){
	dummy = (*pulseIn)[i];
	eIdx2 = i;
      }
    }

    double accumulator = 0;
    for(int i = eIdx1; i < eIdx2 +1; i++){
      accumulator += (*pulseIn)[i]*stepSize;
    }
    
    return accumulator;
    
  }


}




// First we read in the file we're sent from the main function
