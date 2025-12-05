// This is designed for reading in data from our oscilliscope csv
// files for further analysis. Smoothing, derivatives, and
// peak finding algorithms are also included. I want to add a
// trained machine learning algorithm for peak finding to vote on
// how many peaks there are.
// Jonathan Folkerts June 2024

#include<csvParser.hh>

namespace User{
  // returns a vector holding 2d data: Time, ch1, ch2, ch3, ch4
  std::vector<std::vector<double>> tektronicsCSVReader(std::string filename){
    std::vector<std::vector<double>> thisSet;
    std::string line;
    
    std::ifstream file;
    file.open(filename.c_str());
    
    if (!file){
      std::cout << "Error: file \'" + filename + "\' not found." << std::endl;
      std::vector<double> pushme = {-1}; // return -1 so we can error handle
      thisSet.push_back(pushme);
      return thisSet;}
    else{
      for(int j=0;j<21;j++){
	getline(file,line,'\n'); // removes first 21 lines of header
      }
      



      int thisLine = 0;
      int totalColumnNumber = 0;
      
      while(getline(file,line,'\n')){
	// count number of columns
	if(thisLine == 0){
	  std::stringstream sLine;
	  sLine << line.c_str();
	  std::string datum;
	  while(getline(sLine,datum,',')){
	    if(datum != "") totalColumnNumber++;
	  }
	}
	
	std::stringstream sLine;
	sLine << line.c_str();
	std::string datum;
	
       	if(line.empty()){break;}
       	if(line[0] == ',' || line[0] == '\r' || line[0] == '\n' ) {continue;}

	int dataIdx = 0;
	std::vector<double> pushToBack = {};
	while( getline(sLine,datum,',') ){
	  //std::cout << datum << " = datum. Remove after testing.\n";
	  if( datum == "" ) continue;

	  /* left as a comment so I can fix it in the digital version
	  // this needs to be fixed to be more robust
	  if(dataIdx > 4 && AnalogChannelQuantity == 4){
	  dataIdx++;
	  continue;} // digital data has a different time base and breaks things; skip for now
	  */
	  
	  if(datum == "") datum = "-1";// the code doesn't like filler columns, It's dumb and I hate it
	  if(datum == "0.0000e+00") datum = "1e-200";// the code doesn't like 0, It's dumb and I hate it
	  if(datum == "0.000e+00") datum = "1e-200";// the code doesn't like 0, It's dumb and I hate it
	  if(datum == "-0.0000e+00") datum = "1e-200";// the code doesn't like 0, It's dumb and I hate it
	  if(datum == "-0.000e+00") datum = "1e-200";// the code doesn't like 0, It's dumb and I hate it
	  if(datum == "0") datum = "1e-17";// the code doesn't like 0, It's dumb and I hate it
	  if(datum == " -inf") datum = "-INF";// the code doesn't like 0, It's dumb and I hate it
	  if(datum == " inf") datum = "INF";// the code doesn't like 0, It's dumb and I hate it
	  if(datum == "inf") datum = "INF";// the code doesn't like 0, It's dumb and I hate it
	  if(datum == "-inf") datum = "-INF";// the code doesn't like 0, It's dumb and I hate it

	  
	  pushToBack.push_back(std::stod(datum));
	  
	  if(dataIdx < totalColumnNumber){// number of columns to read
	    dataIdx++;
	  }
	  else{
	    dataIdx = 0;
	  }
	}

	thisSet.push_back(pushToBack);

	thisLine ++;
      }
    
    }
    
    file.close();

    // transpose the matrix for use elsewhere
    size_t n = thisSet.size();
    size_t m = thisSet[0].size();
    
    std::vector<std::vector<double>> theReturn(m, std::vector<double>(n));

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            theReturn[j][i] = thisSet[i][j];
        }
    }

    return theReturn;
  }



  // ----------oooooooooo0000000000oooooooooo----------
  
  std::vector<std::vector<double>> tektronicsCSVReaderDigital(std::string filename, int digitalChannelCount){
    std::cout << "std::vector<std::vector<double>> tektronicsCSVReaderDigital(std::string filename, int digitalChannelCount) is not yet implemented. Returning {{-1}}\n";
    std::string line;

    std::vector<std::vector<double>> thisSet;
    std::vector<double> negativeOne = {-1};
    thisSet.push_back(negativeOne);
    return thisSet;
  }



  // ----------oooooooooo0000000000oooooooooo----------


  
  // returns the vertical scales of the csv
  std::vector<double> tektronicsScaleReader(std::string filename){
    std::vector<double> theScales = {}; //empty scales vector
    std::string line;

    std::ifstream file;
    file.open(filename.c_str());
    
    if (!file){
      std::cout << "Error: file \'" + filename + "\' not found." << std::endl;
      std::vector<double> negativeOne = {-1};
      return negativeOne;}
    else{
      for(int j = 0; j < 21; j++){
	if(j != 14 && j != 8){
	  getline(file,line,'\n');}// removes first 21 lines of header
	
	if(j==8){
	  getline(file,line,','); // skip the first column; says "Sample Rate" in csv
	  getline(file,line,','); // Get the double value for the sample rate
	  double pushMe = std::stod(line);
	  theScales.push_back(pushMe);
	  getline(file,line,'\n'); // finish removing line from file
	}

	//Want to read in the vertical scale so we can systematically assign bins
	if( j == 14 ){
	  while( getline(file,line,',') ){
	    //std::cout << line << " = line. Remove after testing.\n";
	    if( line == "Vertical Position" ) continue;
	    if( line == "" ) continue;
	    double pushMe = std::stod(line);
	    theScales.push_back(pushMe);
	    if( line.back() == '\n') break;
	  } // end while
	  //std::cout << "after while\n";
	} // end if
      } // end for
    } // end else

    

    file.close();
    return theScales;
  }
  
  
  
  oscilloscopeData tekToOscopeData(std::string filename) {
    std::ifstream file(filename);
    oscilloscopeData data;
    std::string line;
    int lineCount = 0;
    std::vector<std::vector<double>> channels;
    std::vector<std::vector<double>> digitalChannels;
    bool analog = false, digital = false;

    if (!file.is_open()) {
      std::cerr << "Error opening file: " << filename << std::endl;
      return data;  // Return empty data on failure to open file
    }

    while (getline(file, line)) {
      std::istringstream iss(line);
      std::string key;
      std::getline(iss, key, ',');
      std::string value;
      std::getline(iss, value);

      lineCount++;
      //std::cout << "Line " << lineCount << " = " << line << "\n";
      if (lineCount == 1) {
	data.oscilloscopeModel = value;
      } else if (lineCount == 4) {
	std::istringstream ss(value);
	std::string temp;
	while (getline(ss, temp, ',')) {
	  /*
	  std::cout << temp << '\n';
	  for(char c : temp){
	    std::cout << int(c) << '\n';
	  }
	  */
	  if (temp == "ANALOG" || temp == "ANALOG\r" ) {
	    analog = true;
	  }
	  if (temp == "DIGITAL" || temp == "DIGITAL\r" ) {
	    digital = true;
	    break; // digital is always second
	  }
	}
      } else if (lineCount == 6) {
	if( analog || digital){
	  std::istringstream ss(value);
	  std::string unitStr;

	  // Read and parse the comma-separated values in the line
	  while (getline(ss, unitStr, ',')) {
	    if ( unitStr != "" ) {
	      // The first string found is assumed to be the unit
	      data.horizontalUnit = unitStr;
	      break;
	    }
	  }
	}
      } else if (lineCount == 7) {
	if( analog && digital){
	  std::istringstream ss(value);
	  std::string scaleStr;
	  bool foundAnalogScale = false;

	  // Read and parse the comma-separated values in the line
	  while (getline(ss, scaleStr, ',')) {
            try {
	      double scaleValue = std::stod(scaleStr); // Attempt to convert string to double
	      if (!foundAnalogScale) {
		// If the first double is found, it is assumed to be the analog scale
		data.horizontalScale = scaleValue;
		foundAnalogScale = true;
	      } else {
		// The next double found is assumed to be the digital scale
		data.horizontalDigitalScale = scaleValue;
		break; // Exit the loop after finding the digital scale
	      }
            } catch (const std::invalid_argument& ia) {
	      // If conversion fails, continue to the next CSV
	      continue;
            }
	  }
	}
	else if(analog){
	  data.horizontalScale = std::stod(value);
	}
      } else if (lineCount == 10) {
	//std::cout << "I'm in line 10!\n";
	if (analog || digital) {
	  std::istringstream ss(value);
	  std::string lengthStr;
	  bool foundAnalogLength = false;

	  // Read and parse the comma-separated values in the line
	  while (getline(ss, lengthStr, ',')) {
            try {
	      unsigned long lengthValue = (unsigned long)std::stod(lengthStr); // Attempt to convert string to unsigned long
	      //std::cout << "lengthValue = " << lengthValue << "\n";
	      if (!foundAnalogLength) {
		// If the first valid value is found, it is assumed to be the analog length
		data.recordLengths.push_back(lengthValue);
		foundAnalogLength = true;
	      } else {
		// The next valid value found is assumed to be the digital length
		data.digitalRecordLengths.push_back(lengthValue);
		break; // Exit the loop after finding the digital length
	      }
            } catch (const std::invalid_argument& ia) {
	      // If conversion fails, continue to the next value
	      continue;
            }
	  }
	} else if (analog) {
	  data.recordLengths.push_back(std::stoul(value));
	} else {
	  std::cerr << "No analog data available." << std::endl;
	  break;  // Exit processing as there's no data to process
	}
      } else if (lineCount == 12) {
	std::istringstream ss(value);
	std::string temp = {};
	while (getline(ss, temp, ',') && !temp.empty()) {
	  data.probeAttenuation.push_back(std::stod(temp));
	}
      } else if (lineCount == 13) {
	if (analog && digital) {
	  std::istringstream ss(value);
	  std::string unit;
	  bool switchToDigital = false;

	  // Parse the comma-separated values in the line
	  while (getline(ss, unit, ',')) {
            if (!switchToDigital) {
	      if (unit == "Vertical Units") {
		switchToDigital = true; // Switch to digital units after this point
	      } else {
		data.verticalUnits.push_back(unit);
	      }
            } else {
	      data.digitalVerticalUnits.push_back(unit); // Start filling digital vertical units
            }
	  }
	} else if (analog) {
	  std::istringstream ss(value);
	  std::string unit;
	  while (getline(ss, unit, ',') && !unit.empty()) {
            data.verticalUnits.push_back(unit);
	  }
	} else {
	  std::cerr << "No analog data available." << std::endl;
	  break;  // Exit processing as there's no relevant data to process
	}
      } else if (lineCount == 15) {
	std::istringstream ss(value);
	std::string temp = {};
	while (getline(ss, temp, ',') && !temp.empty()) {
	  data.verticalScale.push_back(std::stod(temp));
	}
      } else if (lineCount == 16) {
	std::istringstream ss(value);
	std::string temp = {};
	while (getline(ss, temp, ',') && !temp.empty()) {
	  data.verticalPosition.push_back(std::stod(temp));
	}
      }  else if (lineCount == 20) {
	std::istringstream ss(value);
	std::string channelLabel = {};
	while (getline(ss, channelLabel, ',') && !channelLabel.empty()) {
	  data.channelLabels.push_back(channelLabel);
	  
	}
      } else if (lineCount == 21) {
	if (analog && digital) {
	  std::istringstream ss(line);
	  std::string channelStr;
	  bool countingDigital = false;
	  int channelCount = 0;
	  while (getline(ss, channelStr, ',')) {
	    if (channelStr == "TIME" || key == "TIME") {
	      key = ""; // get rid of key so we don't get stuck in 
	                // some time nonsense
	      if (countingDigital) {
		data.numberOfDigitalChannels = channelCount;
		break;
	      } else {
		data.numberOfChannels = channelCount;
		countingDigital = true;
		channelCount = 0; // Reset for digital channel counting
	      }
	    } else if (!channelStr.empty()) {
	      channelCount++;
	    } // end if1
	    
	  } // end while
	  channels.resize(data.numberOfChannels);
	  digitalChannels.resize(data.numberOfDigitalChannels);
	} else if (analog){
	  // Determine number of channels by counting commas in a line of data
	  data.numberOfChannels = std::count(line.begin(), line.end(), ',');
	  channels.resize(data.numberOfChannels);
	}
	else{
	  std::cout << "There was not any analog data.\n";
	  break;
	}
	//std::cout << "Check";
      } else if (lineCount > 21) {
	if (lineCount > data.recordLengths[0] + 21) break;
	std::istringstream ss(line);
	std::string timeStr;
	getline(ss, timeStr, ',');
	data.timeData.push_back(std::stod(timeStr));

	for (int i = 0; i < data.numberOfChannels; i++) {
	  std::string channelStr;
	  if (!getline(ss, channelStr, ',')) break;
	  channels[i].push_back(std::stod(channelStr));
	}
      }
    }

    for (auto& channel : channels) {
      data.channelData.push_back(channel);
    }

    file.close();
    return data;
  }



  
  // Smooths the tek data for derivative
  std::vector<double> generateGaussianKernel(int kernelSize, double sigma){
    std::vector<double> kernel(kernelSize);
    double sum = 0.0;
    int halfSize = kernelSize / 2;
    for (int i = 0; i < kernelSize; ++i) {
        double x = i - halfSize;
        kernel[i] = exp(-0.5 * x * x / (sigma * sigma)) / (sqrt(2 * M_PI) * sigma);
	//std::cout << kernel[i] << "\n";
        sum += kernel[i];
    }

    // Normalize the kernel to ensure the sum of all elements equals 1
    for (double& value : kernel) {
        value /= sum;
    }

    return kernel;
  }

  std::vector<double> tekGausSmooth(oscilloscopeData data, size_t channelNumber, double stDev){
    if(stDev == 0){
      return data.channelData[channelNumber];
    }
    double sigmaIndex = 2*stDev; // twice the given standard deviation
    int kernelSize = 3*sigmaIndex; // p/m 3 sigma points around
    int halfSize = kernelSize/2; // auto rounds down; for padding
    std::vector<double> kernel = generateGaussianKernel( kernelSize, sigmaIndex);


    const auto& channel = data.channelData[channelNumber];
    std::vector<double> paddedData(channel.size() + 2 * halfSize);

    // Calculate the average of the first 5 data points
    double averageFirst = 0.0;
    for (int i = 0; i < 5 && i < channel.size(); ++i) {
        averageFirst += channel[i];
    }
    averageFirst /= std::min(5, static_cast<int>(channel.size()));

    // Pad the beginning with the average of the first 5 data points
    std::fill(paddedData.begin(), paddedData.begin() + halfSize, averageFirst);

    // Copy original data to the middle of paddedData
    std::copy(channel.begin(), channel.end(), paddedData.begin() + halfSize);


    // Calculate average of the last 5 data points
    if (channel.size() >= 15) {
      double sumLast15 = 0;
      for (int i = 0; i < 15; ++i) {
        sumLast15 += channel[channel.size() - 15 + i];
      }
      double averageLast15 = sumLast15 / 15;

      // Use the average to fill the extended part of the data
      std::fill(paddedData.begin() + halfSize + channel.size(), paddedData.end(), averageLast15);
    } else {
      // Not enough points to average, so extend with the last point
      std::fill(paddedData.begin() + halfSize + channel.size(), paddedData.end(), channel.back());
    }

    
    std::copy(data.channelData[channelNumber].begin(),
	      data.channelData[channelNumber].end(),
	      paddedData.begin() + halfSize);  // Copy original data to the middle of paddedData

    std::vector<double> smoothedData(data.channelData[channelNumber].size(),
				     0.0);

    for (size_t i = 0; i < smoothedData.size(); ++i) {
      double smoothedValue = 0.0;
      for (int j = -halfSize; j <= halfSize; ++j) {
	int dataIndex = i + j + halfSize;  // Adjust index to access paddedData
	smoothedValue += paddedData[dataIndex] * kernel[j + halfSize];
      }
      smoothedData[i] = smoothedValue;
    }

    return smoothedData;
  }

  

  // Returns the time derivative of the oscilloscope data
  std::vector<double> tekDerivative(oscilloscopeData data, size_t channel){
    std::vector<double> derivative;
    if (channel >= data.numberOfChannels) {
      throw std::out_of_range("Channel index out of range");
    }
    if (data.timeData.size() < 2) {
      throw std::invalid_argument("Not enough data points for differentiation");
    }

    const auto& yValues = data.channelData[channel];
    const auto& xValues = data.timeData;

    derivative.resize(yValues.size(), 0.0);  // Initialize derivative vector with zeros

    
    double h = (xValues.back() - xValues.front()) / (xValues.size() - 1); // time steps are uniform

    // Compute derivatives
    for (size_t i = 0; i < yValues.size(); ++i) {

      if (i == 0) {
	// Forward difference
	derivative[i] = (yValues[i + 1] - yValues[i]) / h;
      } else if (i == yValues.size() - 1) {
	// Backward difference
	derivative[i] = (yValues[i] - yValues[i - 1]) / h;
      } else if (i == 1 || i == yValues.size() - 2) {
	// 2-point central difference
	derivative[i] = (yValues[i + 1] - yValues[i - 1]) / (2 * h);
      } else {
	// 5-point central difference
	derivative[i] = (-yValues[i + 2] + 8 * yValues[i + 1] - 8 * yValues[i - 1] + yValues[i - 2]) / (12 * h);
      }
    }

    return derivative;
  }
  
  
  
  // Returns approximate index of peaks from a smoothed derivatives
  std::vector<std::vector<size_t>> peakLocationsFromDerivative(User::oscilloscopeData derivativeData, User::oscilloscopeData data){
    
    std::vector<std::vector<size_t>> theReturn = {};
    for(size_t i = 0; i < derivativeData.channelData.size(); i++){
      std::vector<size_t> thesePeaks = {};
      // Set absolute threshold for counting a derivative as 8 * the smallest change * timestep
      double derivativeThreshold = 300 * derivativeData.verticalScale[i]/5 / derivativeData.horizontalScale/100;
      // std::cout << " derivativeThreshold = " << derivativeThreshold << "\n";
      bool peakStarting = false;
      size_t peakStart, peakEnd;
      double peakSign = 1;
      bool signSet = false;

      // The sign setting logic can be fooled by the tail of a waveform coming
      // down before the peak from the trigger, but there's not much to do
      // about that. I do pick up on peaks near the tail end which are going
      // up but don't come down.
      
      for(size_t j = 0; j < derivativeData.channelData[i].size(); j++){
	// ignore nans and infs
	if(isnan(derivativeData.channelData[i][j])
	   || isinf(derivativeData.channelData[i][j]) ) {
	  // increase threshold for rebound signal
	  derivativeThreshold = 450 * derivativeData.verticalScale[i]/5 / derivativeData.horizontalScale/100;
	  continue;
	}

	
	// If the data's gotten big enough
	if( !peakStarting && (peakSign*abs(derivativeData.channelData[i][j]) > derivativeThreshold) ){
	  // Check to see if we need to set the sign
	  if(!signSet){
	    signSet = true;
	    peakSign = derivativeData.channelData[i][j]/abs(derivativeData.channelData[i][j]);
	    //std::cout << "Peak Sign is : " << peakSign << "\n";
	  }
	  // skip if original data on the wrong sign
	  if(peakSign * data.channelData[i][j] < 0) continue;
	  
	  // skip wrong sign
	  if(peakSign * derivativeData.channelData[i][j] < 0) continue;
	  //std::cout << "\nPeak is starting at index : " << j << "\n";
	  peakStart = j;
	  peakStarting = true;	  continue; // don't need the rest of the loop
	}

	
	// If we're crossing zero after going up
	if(peakStarting && (peakSign*derivativeData.channelData[i][j] < -derivativeThreshold) ){
	  //std::cout << "Peak is ending at index : " << j << "\n";
	  peakEnd = j;
	  peakStarting = false;
	  //std::cout << "Found a peak at idx " << j << " Where the derivative is: " <<  derivativeData.channelData[i][j] << "\n";
	  thesePeaks.push_back( (peakStart+peakEnd) / 2 );
	  continue; // don't need the rest of the loop
	}
      }
      
    
      // If we had id'd the start of a peak, the end must be somewhere further on
      // not working right currently. It  gets fooled too easily by slight upturns at the end from randomness
      if(peakStarting){
	thesePeaks.push_back(derivativeData.channelData.size() + 1);
      }
      theReturn.push_back(thesePeaks);
    }

    //std::cout << "There were " << theReturn[0].size() << " peaks.\n";
    return theReturn;
  }


  
  bool peakPassesEnergyCutIntegral(size_t peakIdx, double min_in_Vs, double max_in_Vs, User::oscilloscopeData data, size_t channel){
    size_t startIdx = peakIdx;
    while(data.channelData[channel][peakIdx] > 5*data.verticalMinorTic[channel]){
      if(startIdx == 0){
	std::cout << "Reached a start index of 0; this peak will have too low energy\n";
	break;
      }
      startIdx--;
    }

    size_t endIdx = peakIdx;
    while(data.channelData[channel][peakIdx] > 5*data.verticalMinorTic[channel]){
      if(startIdx == data.channelData[channel].size()-1){
	std::cout << "Reached an end index of the last entry; this peak will have too low energy\n";
	break;
      }
      endIdx++;
    }

    double integral = 0;
    double timeStep = data.timeData[1] - data.timeData[0];
    for(size_t i = startIdx; i <= endIdx; i++){
      integral += timeStep * data.channelData[channel][i];
    }
    
    bool belowMax = integral < max_in_Vs;
    bool aboveMin = integral > min_in_Vs;

    return aboveMin && belowMax;
  }

  bool peakPassesEnergyCutVoltage(size_t peakIdx, double min_in_V, double max_in_V, User::oscilloscopeData data, size_t channel){
    bool belowMax = data.channelData[channel][peakIdx] < max_in_V;
    bool aboveMin = data.channelData[channel][peakIdx] > min_in_V;
    return aboveMin && belowMax;
  }
  
  
}

