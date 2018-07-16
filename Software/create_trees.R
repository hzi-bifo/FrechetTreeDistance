#!/usr/bin/env Rscript

# This script takes a multiple sequence alignment and a list of locations for each sequence as input.
# It infers five different phylogentic trees using different methods (Parsimony (Fitch algorithm), NJ, UPGMA, ML (JC and GTR substitution model))
# and does a phylogeographic reconstruction using parsimony.

# Run this script from the command line using:
# Rscript create_trees.R [alignment] [locations] [outgroup for rooting]

### define functions ###

#ancestral character state reconstruction using parsimony
ancestral_reconstruction <- function(phylo, locdat){
  ancestors_loc <- ancestral.pars(phylo, locdat, type = "ACCTRAN") #without cost matrix specified: Fitch algorithm
  
  n_tips <- Ntip(phylo) #number of tips
  n_nodes <- length(unique(phylo$edge[,1])) #number of internal nodes
  n <- n_tips + n_nodes
  
  all_loc <- levels(locdat)
  node_locs <- vector(mode="character", length=n)
  
  for (i in 1:n) {
    node_locs[i] <- all_loc[which(ancestors_loc[[i]] == max(ancestors_loc[[i]]))][1]
  }
  
  # add internal node IDs
  node_label <- paste("intNode", seq(1, n_nodes), sep="")
  phylo$node.label <- node_label
  
  # create annotation
  node_annotation <- data.frame(label=c(phylo$tip.label, phylo$node.label), location=node_locs)
  
  return(list(tree=phylo, annotation=node_annotation))
}

infer_tree <- function(seq_data, method, model="none"){
  # get distances
  dist <- dist.ml(seq_data) #standard model: JC69, alternative is F81
  
  if (method == "UPGMA"){
    phylo <- upgma(dist) #rooted
  } else if (method == "NJ"){
    phylo <- NJ(dist) #unrooted
  } else if (method == "Parsimony"){
    start <- upgma(dist) #rooted, UPGMA as starting tree
    phylo <- optim.parsimony(start, seq_data) #optimize using parsimony criterium
  } else if (method == "ML"){
    start <- upgma(dist) #rooted, UPGMA as starting tree
    ml <- pml(start, data=seq_data) #get likelihood
    if (model == "JC"){
      ml_result <- optim.pml(ml, TRUE) #standard parameters: Jukes Cantor model
    } else if(model == "GTR"){
      ml_result <- optim.pml(ml, TRUE, model="GTR")
    } else {
      ## throw error
    }
    phylo <- ml_result$tree #unrooted
  } else {
    ## throw error
  }
  return(phylo)
}

root_tree <- function(phylo, outgroup){
  # roots tree using an outgroup and subsequently removes outgroup
  phylo_rooted <- root(phylo, outgroup = outgroup, resolve.root=TRUE)
  phylo_rmoutgroup <- drop.tip(phylo_rooted, outgroup)
  return(phylo_rmoutgroup)
}

save_tree <- function(result, filename){
  write.tree(result$tree, paste(filename,".phy", sep=""))
  write.table(result$annotation, paste(filename, ".annotation.txt", sep=""), quote=FALSE, row.names=FALSE, sep="\t")
}

### perform analysis ###

library("phangorn") # used for tree inference and ancestral character state reconstruction

args = commandArgs(trailingOnly=TRUE)
alignment_file <- args[1]
locations <- args[2]
outgroup <- args[3]

# load alignment as phyDat object
seq_data <- read.phyDat(alignment_file, format = "fasta")

# tree inference: parsimony, distance matrix, maximum likelihood
tree_UPGMA <- infer_tree(seq_data, "UPGMA")
tree_NJ <- infer_tree(seq_data, "NJ")
tree_Fitch <- infer_tree(seq_data, "Parsimony")
tree_MLJC <- infer_tree(seq_data, "ML", "JC")
tree_MLGTR <- infer_tree(seq_data, "ML", "GTR")

# root trees and remove outgroup
rooted_UPGMA <- root_tree(tree_UPGMA, "f0dp0")
rooted_NJ <- root_tree(tree_NJ, "f0dp0")
rooted_Fitch <- root_tree(tree_Fitch, "f0dp0")
rooted_MLJC <- root_tree(tree_MLJC, "f0dp0")
rooted_MLGTR <- root_tree(tree_MLGTR, "f0dp0")

tree_list <- c(rooted_UPGMA, rooted_NJ, rooted_Fitch, rooted_MLJC, rooted_MLGTR)
names <- c("UPGMA", "NJ", "Parsimony", "MLJC", "MLGTR")
n_trees <- length(tree_list)

# for information: check if trees have the same topology

test.topology <- matrix(data = NA, nrow=n_trees, ncol=n_trees)
rownames(test.topology) <- names
colnames(test.topology) <- names

for (i in 1:n_trees){
  for (j in 1:n_trees){
    #calculate distances using tree topology
    test.topology[i,j] <- RF.dist(tree_list[[i]], tree_list[[j]])
  }
}

write.table(test.topology, "pairwiseRFDistances.txt", quote=FALSE, sep="\t")

# add ancestral character states

location_info <- read.table(locations, header=TRUE, sep="\t")

#observed locations
tipdata <- as.matrix(location_info$location)
row.names(tipdata) <- location_info$id
 
#save observed locations in phyDat format
all_locs <- unique(tipdata)
locdat <- phyDat(tipdata, type="USER", levels=all_locs)

# ancestral state reconstruction
dir.create(file.path(getwd(),"trees"), showWarnings = FALSE) # creates output directory if it doesn't exist

ancestral_UPGMA <- ancestral_reconstruction(rooted_UPGMA, locdat)
save_tree(ancestral_UPGMA, "trees/UPGMA")

ancestral_NJ <- ancestral_reconstruction(rooted_NJ, locdat)
save_tree(ancestral_NJ, "trees/NJ")

ancestral_Fitch <- ancestral_reconstruction(rooted_Fitch, locdat)
save_tree(ancestral_Fitch, "trees/Parsimony")

ancestral_MLJC <- ancestral_reconstruction(rooted_MLJC, locdat)
save_tree(ancestral_MLJC, "trees/MLJC")

ancestral_MLGTR <- ancestral_reconstruction(rooted_MLGTR, locdat)
save_tree(ancestral_MLGTR, "trees/MLGTR")