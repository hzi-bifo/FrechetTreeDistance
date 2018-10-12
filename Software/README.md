## Running the pipeline

The script [Pipeline.sh](https://github.com/hzi-bifo/FrechetTreeDistance/blob/master/Software/Pipeline.sh) performs the tree inference and phylogeographic reconstruction, calculates the discrete Fréchet distances and visualizes them using multidimensional scaling.

### 1) Requirements

The pipeline is implemented in bash and combines different R scripts. Necessary software packages (and versions this software was tested for) are:
* R and Rscript (3.4.4)
* R packages: ape (5.1), phangorn (2.4.0), ggplot2 (2.2.1), MASS (7.3-50)

### 2) Input

The input files used in this manuscript are located in the [Data/Raw directory](https://github.com/hzi-bifo/FrechetTreeDistance/tree/master/Data/Raw).

The necessary input for the pipeline is:
* [A multiple sequence alignment](https://github.com/hzi-bifo/FrechetTreeDistance/blob/master/Data/Raw/Wallace_H5N1_HA.aln) in fasta format. It needs to be named *\[PREFIX\].aln*.
* [A tab delimited file](https://github.com/hzi-bifo/FrechetTreeDistance/blob/master/Data/Raw/Wallace_H5N1_HA.locations.txt) with one column called "id" and one column called "location" listing the sequence identifier and the assigned location for the sequence. It needs to be named *\[PREFIX\].locations.txt*.
* [A csv file](https://github.com/hzi-bifo/FrechetTreeDistance/blob/master/Data/Raw/distance.matrix.csv) containing the symmetric distance matrix with distances between all observed locations and including the locations as row and column names.
* The identifier of the sequence that that will be used to root all trees. This sequence should be an outgroup and is removed from the dataset after rooting.

To run the pipeline on different input data, edit the lines 5 (prefix), 6 (distance matrix) and 7 (root sequence) in [Pipeline.sh](https://github.com/hzi-bifo/FrechetTreeDistance/blob/master/Software/Pipeline.sh).

### 3) Running the pipeline 

To run the analysis and replicate all results from the manuscript, change into the folder containing the [raw data](https://github.com/hzi-bifo/FrechetTreeDistance/tree/master/Data/Raw) and run
> bash path/to/folder/Pipeline.sh

### 4) Output

The pipeline produces the following output files in the directory where the script was called:
* Five phylogenetic trees inferred using different methods (*\[method\].phy*) in Newick format in a folder called trees. All nodes in these trees are labeled.
* For each tree the corresponding phylogeographic reconstruction (*\[method\].annotation.txt*) given in a tab delimited file with a column called "label" containing the node label and a column called "location" containing the inferred location (in case of internal nodes) or the given location (in case of leaf nodes).
* A text file (*pairwiseRFDistances.txt*) containing Robinson-Foulds distances between all trees.
* A text file (*pairwiseFrechetDistances.txt*) containing discrete Fréchet tree distances between all trees.
* A png image (*pairwiseFrechetDistances.png*) containing the multidimensional scaling plot.

## Calculating the Fréchet tree distance

Instead of running the complete pipeline, the discrete Fréchet tree distance can be calculated on arbitrary input trees. The distance measure is implemented in the script [FrechetTreeDistance.R](https://github.com/hzi-bifo/FrechetTreeDistance/blob/master/Software/FrechetTreeDistance.R).

### 1) Requirements

Necessary software packages (and versions this software was tested for) are:
* R and Rscript (3.4.4)
* R packages: ape (5.1)

### 2) Input

The necessary input for this script is:
* A phylogenetic tree (named *\[prefix1\].phy*) in Newick format with labels for all nodes.
* A phylogeographic reconstruction for this phylogenetic tree (named *\[prefix1\].annotation.txt*) given in a tab delimited file with a column called "label" containing the node label and a column called "location" containing the inferred location (in case of internal nodes) or the given location (in case of leaf nodes).
* A second phylogenetic tree (named *\[prefix2\].phy*) in Newick format with labels for all nodes.
* A phylogeographic reconstruction for this second phylogenetic tree (named *\[prefix2\].annotation.txt*) given in a tab delimited file with a column called "label" containing the node label and a column called "location" containing the inferred location (in case of internal nodes) or the given location (in case of leaf nodes).
* A csv file containing the symmetric distance matrix with distances between all observed locations and including the locations as row and column names.

The two trees that are compared should be inferred for the same taxa. Both the number of leaves and the label at the leaves need to be equal.

Further optional arguments are possible:
* Adjustment: defines whether distances should be adjusted by the number of paths, i.e. if the cost for each node should be divided by the number of descendant leaves (see the manuscript for details). Needs to be TRUE or FALSE. If the argument is not defined, the default is TRUE. The adjustment is recommended to aid interpretability.

### 3) Running the script

In the folder where the data is located, the script can be called from the command line using:
> Rscript path/to/folder/FrechetTreeDistance.R \[prefix1\] \[prefix2\] \[distance matrix\] \[adjustment\]

The discrete Fréchet tree distance between the two given trees is output on the command line.
