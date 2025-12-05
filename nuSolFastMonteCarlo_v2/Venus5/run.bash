#!/bin/bash

# Ensure the data directory exists
mkdir -p data

# Count the number of .root files already in the data directory
start_index=$(ls data/*.root 2> /dev/null | wc -l)

# Calculate the ending index as start_index + 1000
end_index=$((start_index + 1000))

# Loop from the starting index to the end index
for ((i=start_index; i<end_index; i++))
do
  echo "Running simulation for run ${i}"
  # Execute the Monte Carlo simulation
  ./nuSolFastMonteCarlo -b batch.mac

  # It assumes that each execution generates exactly one .root file in the current directory
  # Move and rename the .root files
  mv *.root "data/run${i}.root"
done

echo "All tasks completed."
