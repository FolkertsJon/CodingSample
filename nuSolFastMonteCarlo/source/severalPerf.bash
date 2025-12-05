#!/bin/bash
# First tests of a shell script to do a bunch of runs at once


for i in 1 2 3 4 5 6 7 8 9 10
do
    name = 'output"$i".root'
    echo 'gROOT->LoadMacro("nuSolPerformance.cpp"); gROOT->Macro("nuSolPerformance.cpp"); gSystem->Exit(0);' | root -b -l
    
    mv outfile.root 'output'"$i"'.root'
done
