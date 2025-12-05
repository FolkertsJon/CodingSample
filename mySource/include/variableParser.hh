// This is the header file for variableParcer.cpp
// Jonathan Folkerts May 2021
#ifndef variableParser_H
#define variableParser_H

#include<algorithm>
#include<fstream>
#include<string>
#include<sstream>
#include<iostream>
#include<vector>
#include<cstring>

#include<kUnits.hh>
#include<galliumInteraction.hh>
#include<backgrounds.hh>
#include<dataStructures.hh>

namespace User
{
  extern bool printedHelp;
  void parInitializer(User::Par *p);// Fills default paramters
  void parChanger(User::Par *p, User::Variable v);
  char* getCmdOption(char ** begin, char ** end, const std::string & option);
  bool cmdOptionExists(char** begin, char** end, const std::string& option);
  bool isComment(std::string theString);
  void variablePusher(std::vector<User::Variable> &holder, std::string file);
  void printHelp();
  int argumentInterpreter(int argc, char * argv[], User::Par *par);// Takes arguments from main and uses them
}

#endif
