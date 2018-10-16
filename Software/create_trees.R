#!/usr/bin/env Rscript

# This script takes a multiple sequence alignment and a list of locations for each sequence as input.
# It infers five different phylogentic trees using different methods (Parsimony (Fitch algorithm), NJ, UPGMA, ML (JC and GTR substitution model))
# and does a phylogeographic reconstruction.

# Run this script from the command line using:
# Rscript create_trees.R [alignment] [locations] [outgroup for rooting] [method for ancestral state reconstruction]
# The last argument is optional. Default is "Parsimony". Possible values are "Parsimony" and "ML".

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

ml_ancestral_reconstruction <- function(phylo, tipdata){
  x <- tipdata[phylo$tip.label,]
  
  #function doesn't allow values of 0
  phylo$edge.length[phylo$edge.length <= 0] <- 1e-100
  ancestral_ace <- ace(x, phylo, type="discrete", method = "ML")
  #ancestral_ace <- reconstruct(x, phylo, method = "ML")
  
  likelihoods <- ancestral_ace$lik.anc
  n_nodes <- nrow(likelihoods)
  node_locs <- vector(mode="character", length=n_nodes)
  for (i in 1:n_nodes){
    node_locs[i] <- names(which(likelihoods[i,] == max(likelihoods[i,])))
  }
  
  # add internal node IDs
  node_label <- paste("intNode", seq(1, n_nodes), sep="")
  phylo$node.label <- node_label
  
  # create annotation
  node_annotation <- data.frame(label=c(phylo$tip.label, phylo$node.label), location=c(x, node_locs))
  
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
    phylo <- acctran(phylo, seq_data) #get branch lengths using acctran
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
  # write.tree replaces spaces in sequence names by underscores - do the same with annotation file
  result$annotation$label <- gsub(" ", "_", result$annotation$label)
  write.table(result$annotation, paste(filename, ".annotation.txt", sep=""), quote=FALSE, row.names=FALSE, sep="\t")
}

### perform analysis ###

library("phangorn") # used for tree inference and ancestral character state reconstruction

args = commandArgs(trailingOnly=TRUE)
alignment_file <- args[1]
locations <- args[2]
outgroup <- args[3]
asr_method <- args[4]

# default for ancestral state reconstruction: Parsimony
if (is.na(asr_method)){
  asr_method = "Parsimony"
}
if (asr_method != "Parsimony" && asr_method != "ML"){
  stop("Invalid argument. Method for ancestral state reconstruction must be set to 'Parsimony' or 'ML'")
}

# load alignment as phyDat object
seq_data <- read.phyDat(alignment_file, format = "fasta")

# test if outgroup is present in alignment
if (!outgroup %in% names(seq_data)){
  stop(paste('The outgroup', outgroup, 'is not present in the aligned sequence file.'))
}

# tree inference: parsimony, distance matrix, maximum likelihood
tree_UPGMA <- infer_tree(seq_data, "UPGMA")
tree_NJ <- infer_tree(seq_data, "NJ")
tree_Fitch <- infer_tree(seq_data, "Parsimony")
tree_MLJC <- infer_tree(seq_data, "ML", "JC")
tree_MLGTR <- infer_tree(seq_data, "ML", "GTR")

# root trees and remove outgroup
rooted_UPGMA <- root_tree(tree_UPGMA, outgroup)
rooted_NJ <- root_tree(tree_NJ, outgroup)
rooted_Fitch <- root_tree(tree_Fitch, outgroup)
rooted_MLJC <- root_tree(tree_MLJC, outgroup)
rooted_MLGTR <- root_tree(tree_MLGTR, outgroup)

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

# test if ids in location info match ids in sequence file
for (seq in rooted_UPGMA$tip.label){
  if (! seq %in% location_info$id){
    stop(paste('The sequence', seq, 'is not present in the location file.'))
  }
}

#observed locations
tipdata <- as.matrix(location_info$location)
row.names(tipdata) <- location_info$id
 
#save observed locations in phyDat format
all_locs <- unique(tipdata)
locdat <- phyDat(tipdata, type="USER", levels=all_locs)

# ancestral state reconstruction
dir.create(file.path(getwd(),"trees"), showWarnings = FALSE) # creates output directory if it doesn't exist

if (asr_method == "Parsimony"){
  ancestral_UPGMA <- ancestral_reconstruction(rooted_UPGMA, locdat)
  ancestral_NJ <- ancestral_reconstruction(rooted_NJ, locdat)
  ancestral_Fitch <- ancestral_reconstruction(rooted_Fitch, locdat)
  ancestral_MLJC <- ancestral_reconstruction(rooted_MLJC, locdat)
  ancestral_MLGTR <- ancestral_reconstruction(rooted_MLGTR, locdat)
} else if (asr_method == "ML"){
  ancestral_UPGMA <- ml_ancestral_reconstruction(rooted_UPGMA, tipdata)
  ancestral_NJ <- ml_ancestral_reconstruction(rooted_NJ, tipdata)
  ancestral_Fitch <- ml_ancestral_reconstruction(rooted_Fitch, tipdata)
  ancestral_MLJC <- ml_ancestral_reconstruction(rooted_MLJC, tipdata)
  ancestral_MLGTR <- ml_ancestral_reconstruction(rooted_MLGTR, tipdata)
}

save_tree(ancestral_UPGMA, "trees/UPGMA")
save_tree(ancestral_NJ, "trees/NJ")
save_tree(ancestral_Fitch, "trees/Parsimony")
save_tree(ancestral_MLJC, "trees/MLJC")
save_tree(ancestral_MLGTR, "trees/MLGTR")
