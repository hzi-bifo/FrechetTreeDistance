#!/bin/bash

# set input files and parameters

name="Wallace_H5N1_HA"
distances="distance.matrix.csv"
outgroup="f0dp0"

software_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

# set name for output file and delete if it already exists (results are appended)
distancefile="pairwiseFrechetDistances.txt"
rm -f $distancefile

# create trees and perform ancestral character state reconstruction
# creates .phy and .annotation.txt files for Parsimony, UPGMA, NJ, MLJC, MLGTR
echo "----- 1 Create trees and infer ancestral character states -----"
echo
Rscript $software_dir"/create_trees.R" $name".aln" $name".locations.txt" $outgroup

# calculate pairwise frechet distances
# the script gets called for each pairwise comparison and results are written to the command line as well as the output file
echo " ----- 2 Calculate Fréchet distances -----"

echo -e "Parsimony\tUPGMA\t$(Rscript $software_dir"/FrechetTreeDistance.R" "trees/Parsimony" "trees/UPGMA" $distances)" | tee -a $distancefile
echo -e "Parsimony\tNJ\t$(Rscript $software_dir"/FrechetTreeDistance.R" "trees/Parsimony" "trees/NJ" $distances)" | tee -a $distancefile
echo -e "Parsimony\tMLJC\t$(Rscript $software_dir"/FrechetTreeDistance.R" "trees/Parsimony" "trees/MLJC" $distances)" | tee -a $distancefile
echo -e "Parsimony\tMLGTR\t$(Rscript $software_dir"/FrechetTreeDistance.R" "trees/Parsimony" "trees/MLGTR" $distances)" | tee -a $distancefile

echo -e "UPGMA\tNJ\t$(Rscript $software_dir"/FrechetTreeDistance.R" "trees/UPGMA" "trees/NJ" $distances)" | tee -a $distancefile
echo -e "UPGMA\tMLJC\t$(Rscript $software_dir"/FrechetTreeDistance.R" "trees/UPGMA" "trees/MLJC" $distances)" | tee -a $distancefile
echo -e "UPGMA\tMLGTR\t$(Rscript $software_dir"/FrechetTreeDistance.R" "trees/UPGMA" "trees/MLGTR" $distances)" | tee -a $distancefile

echo -e "NJ\tMLJC\t$(Rscript $software_dir"/FrechetTreeDistance.R" "trees/NJ" "trees/MLJC" $distances)" | tee -a $distancefile
echo -e "NJ\tMLGTR\t$(Rscript $software_dir"/FrechetTreeDistance.R" "trees/NJ" "trees/MLGTR" $distances)" | tee -a $distancefile

echo -e "MLJC\tMLGTR\t$(Rscript $software_dir"/FrechetTreeDistance.R" "trees/MLJC" "trees/MLGTR" $distances)" | tee -a $distancefile

# create multidimensional scaling to plot Frechet distances
Rscript $software_dir"/visualize_distances.R" $distancefile
