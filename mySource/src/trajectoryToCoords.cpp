
#include "trajectoryToCoords.hh"

namespace User
{
double* trajectoryTot(std::string theFile){
  std::ifstream data;
  std::ofstream outData;
  std::string line;
  int nLines = number_of_lines(theFile);
  double datum;
  double *t = new double[nLines];


  data.open(theFile);
  outData.open("fileOut.txt");
  std::cout << "File is open m8.\n";
  int currentLine = 0;
  while(data >> line){
    //std::cout << "The line is: \"" << line << "\"\n";
    std::stringstream s_stream(line);
    bool first = true;
    //std::cout << "The boolean value of first is: " << first << std::endl;
    while( s_stream.good() ){
      std::string substr;
      getline(s_stream, substr, ',');
      //std::cout << "Code says the substring is: \"" << substr << "\"\n";
      if(first){
	std::string::size_type *ptr = 0;
	//std::cout << "Before t the string is " << substr <<" \n";
	double placeholder = stod(substr,ptr);
	t[currentLine] = placeholder;
	//std::cout << "After t\n;";
	first = false;
      }

    }
    currentLine++;
    
    //std::cout << (line) << std::endl;
    outData << (line) << std::endl;
  }
  for(int i=0 ; i<5;i++){
    std::cout << "t[" << i << "]=" << t[i] << ".\n" ;
  }

  data.close();
  outData.close();
  return t;
}

double* trajectoryToR(std::string theFile){
  std::ifstream data;
  std::ofstream outData;
  std::string line;
  int nLines = number_of_lines(theFile);
  double datum;
  double *R = new double[nLines];


  data.open(theFile);
  outData.open("fileOut.txt");
  std::cout << "File is open m8.\n";
  int currentLine = 0;
  while(data >> line){
    //std::cout << "The line is: \"" << line << "\"\n";
    std::stringstream s_stream(line);
    bool first = true;
    //std::cout << "The boolean value of first is: " << first << std::endl;
    while( s_stream.good() ){
      std::string substr;
      getline(s_stream, substr, ',');
      //std::cout << "Code says the substring is: \"" << substr << "\"\n";
      if(first){
	first = false;
      }
      else{
	std::string::size_type *ptr = 0;
	//std::cout << "Before R\n";
	double placeholder = stod(substr,ptr);
	R[currentLine] = 215*placeholder;
	//std::cout << "After R\n";
	first = true;
      }

    }
    std::string::size_type *ptr = 0;
    //R[currentLine] = 215*stod(line,ptr);
    currentLine++;
    
    //std::cout << (line) << std::endl;
    outData << (line) << std::endl;
  }
  for(int i=0 ; i<5;i++){
    std::cout << "R[" << i << "]=" << R[i] << ".\n" ;
  }

  data.close();
  outData.close();
  return R;
}
}
