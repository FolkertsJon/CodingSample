#!/bin/bash

closest_values=(3 5 10 15)  # Replace this with the values you want to iterate over
furthest_values=(81.63 154.8)  # Replace this with the furthest values
#closest_values=(3)  # Replace this with the values you want to iterate over
#furthest_values=(81.63)  # Replace this with the furthest values
modes=(1 0)  # modes (0 =  elliptical; 1 = file)

checkedFile=0

for mode in "${modes[@]}"
do
    for closest in "${closest_values[@]}"
    do
	for furthest in "${furthest_values[@]}"
	do
	    
	    if [ "$checkedFile" = "1" ]; then
		break
	    elif [ "$mode" = "1" ]; then
		checkedFile=1
	    fi
	    
            #Write the string into the "batch.mac" file
            echo "##################################################
# This is the batch file for nuSolPerformance
# a \"#\" deonotes the line is a comment/unused
# lines cannot be empty
##################################################
#
# File (1) or Elliptical (0) mode
fileMode $modes
#
#
# time step in sec
timeStep 0.0001
# 86397 - slightly less than a day so that we get close for year-long runs (at earth)
#
# Detector Radius in Meters
detectorRadius 0.3
#
#
# gallium Mass
kgGallium 100
# 173.76 - 400 kg GaGG + 400 kg Tungsten-Gallium-Phosphate
# 90.27 - 400 kg GaGG 
# For ellipses
# Closest approach in solar radii
closest $closest
# Furthest in solar radii (154.8 = Venus; 81.63 = Mercury)
furthest $furthest
#
#
#
# Number of flights to run
nLoops 1" > batch.mac

            # Execute the program "./nuSolPerformance -b batch.mac"
            ./nuSolPerformance -b batch.mac

            # Check if furthest is 154.8
            if [ "$furthest" = "154.8" ]; then
		# Move "radiusHist.png" to "Venus-closest.png" with the value of closest
		mv radiusHist.png "Venus-$closest.png"
            elif [ "$furthest" = "81.63" ]; then
		# Move "radiusHist.png" to "Mercury-closest.png" with the value of closest
		mv radiusHist.png "Mercury-$closest.png"
            elif [ "$mode" = "1" ]; then
		# Move "radiusHist.png" to "Mercury-closest.png" with the value of closest
		mv radiusHist.png "File.png"
            fi
	done
    done
done

