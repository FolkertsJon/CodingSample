
#include<fstream>
#include<string>
#include<iostream>
#include<sstream>

using namespace std;
  // get number of lines  
  int number_of_lines (std::string theFile) {
    // char* intermediate = theFile;
    std::ifstream intermediate (theFile);
    //data.open (intermediate);
    std::cout << "I have opened the file \"" << theFile << "\"\n\n";
    std::string line;
    int nLines=0;
    while (std::getline(intermediate, line))
      ++nLines;
    std::cout << "Number of lines in text file: " << nLines << std::endl;
    return nLines;
  }


int main(){
  string theFile = "ACOorbits.txt";
  ifstream data;
  ofstream outData;
  string line;
  int nLines = number_of_lines(theFile);
  double datum;
  double R[nLines];
  double t[nLines];


  data.open(theFile);
  outData.open("play.txt");
  cout << "File is open m8.\n";
  int currentLine = 0;
  while(data >> line){
    cout << "The line is: \"" << line << "\"\n";
    stringstream s_stream(line);
    bool first = true;
    cout << "The boolean value of first is: " << first << endl;
    while( s_stream.good() ){
      std::string substr;
      getline(s_stream, substr, ',');
      cout << "Code says the substring is: \"" << substr << "\"\n";
      if(first){
	std::string::size_type *ptr = 0;
	cout << "Before t\n";
	t[currentLine] = stod(substr,ptr);
	cout << "After t\n";
	first = false;
      }
      else{
	std::string::size_type *ptr = 0;
	cout << "Before x\n";
	double placeholder = stod(substr,ptr);
	R[currentLine] = placeholder;
	cout << "After x\n";
	first = true;
      }

    }
    currentLine++;
    double R = 0;
    
    cout << (line) << endl;
    outData << (line) << endl;
  }
  for(int i=0 ; i<100;i++){
    cout << "R[" << i << "]=" << R[i] << " and " 
	 << "t[" << i << "]=" << t[i] << ".\n" ;
  }

  data.close();
  outData.close();


  return 0;
}
