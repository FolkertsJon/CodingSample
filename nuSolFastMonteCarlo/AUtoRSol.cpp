#include<fstream>
#include<string>
#include<iostream>

using namespace std;


int main(){
  ifstream data;
  ofstream outData;
  double datum;
  
  data.open("WSUorbits.txt");
  outData.open("WSUOrbitsRSol.txt");
  cout << "File is open m8.\n";
  while(data >> datum){
    outData << (215*datum) << endl;
  }
  data.close();
  outData.close();


  return 0;
}
