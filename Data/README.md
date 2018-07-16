This folder contains both raw data used to start the analysis and result files.

## Raw data

These files include:
* Wallace_H5N1_HA.fa: fasta file containing the original unaligned sequence data as downloaded from the NCBI influenza database.
* Wallace_H5N1_HA.aln: fasta file containing the aligned and curated sequence data. To avoid potential problems with spaces in the fasta header, all header lines were replaced by a unique identifier.
* Wallace_H5N1_HA.map: text file containing information which new identifier represents which sequence.
* Wallace_H5N1_HA.locations.txt: text file with the assigned location for each sequence identifier.
* distance.matrix.csv: csv file containing the symmetric distance matrix with distances between all observed locations and including the locations as row and column names.

## Result data

These files include:
* Five phylogenetic trees ([method].phy) in Newick format in a folder called trees. All nodes in these trees are labeled.
* For each tree the corresponding phylogeographic reconstruction ([method].annotation.txt) given in a tab delimited file with a column called "label" containing the node label and a column called "location" containing the inferred location (in case of internal nodes) or the given location (in case of leaf nodes).

Additional results like tables with the Robinson Foulds and Fréchet tree distances and figures of all phylogenetic trees can be found in the folder [Figures&Tables](https://github.com/hzi-bifo/FrechetTreeDistance/tree/master/Figures%26Tables).
