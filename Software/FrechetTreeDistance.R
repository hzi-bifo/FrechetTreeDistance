#!/usr/bin/env Rscript

# This script implements the discrete Frechet tree distance.

# Run this script from the command line using:
# Rscript FrechetTreeDistance.R [name of reference] [name of test] [distance matrix]

### define functions ###

read_tree <- function(filename) {
  # read tree and annotation (format as in output from create_trees.R)
  tree <- read.tree(paste(filename, ".phy", sep=""))
  annotation <- read.table(paste(filename, ".annotation.txt", sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE)
  return(list(tree, annotation))
}

# takes a tree and a tip (represented by the label) and returns the path (both locations and node ids) from the root to this node
get_path <- function(tree, tip, tip_locations) {
  # get start id
  start <- which(tree$tip.label == tip);
  # get start location
  start_loc <- tip_locations$location[which(tip_locations$label == tip)]
  
  # initialize character vector to save path of locations
  path <- character();
  # initialize integer vector to save ids of nodes in this path
  path_id <- integer();
  
  # save all locations in one vector
  all_label <- c(tree$tip.label, tree$node.label)
  
  # start at tip and go down the tree until the root node is reached
  stop <- 0;
  id <- start;
  
  # set first location
  path <- c(start_loc, path);
  path_id <- c(id, path_id);
  
  while(stop == 0) {
    # look for incoming edges
    position <- which(tree$edge[,2] == id);
    # if no incoming edge: current node is the root (leave the while loop)
    if (length(position) == 0) {
      stop <- 1;
    }
    # if there is an incoming edge: get id of next node
    id <- tree$edge[position, 1];
    # add location and id
    path <- c(tip_locations$location[which(tip_locations$label == all_label[id])], path);
    path_id <- c(id, path_id);
  }
  return(list("locations" = path, "ids" = path_id));
}

# computes Frechet coupling to calculate the distance between two paths of locations in the next step
frechetCoupling <- function(path1, path2, distance_matrix, sum_or_max = "sum") {
  p <- length(path1)
  q <- length(path2)
  
  # initialize empty matrices
  distances <- matrix(0, p, q)
  frechet <- matrix(0, p, q)
  
  # fill frechet matrix
  for (i in 1:p) {
    for (j in 1:q) {
      # get distance between the current locations
      distances[i, j] <- distance_matrix[path1[i], path2[j]]
      if (i == 1 && j == 1) {
        frechet[1, 1] = distances[1, 1]
      }
      if (i > 1 && j == 1) {
        frechet[i, 1] = do.call(sum_or_max, list(frechet[i - 1, 1], distances[i, 1]))
      }
      if (i == 1 && j > 1) {
        frechet[1, j] = do.call(sum_or_max, list(frechet[1, j - 1], distances[1, j]))
      }
      if (i > 1 && j > 1) {
        frechet[i, j] = do.call(sum_or_max, list(min(frechet[i - 1, j], frechet[i - 1, j - 1], frechet[i, j - 1]), distances[i, j]))
      }
    }
  }
  
  # traceback to determine which points are connected
  seq1 <- i <- p
  seq2 <- j <- q
  
  while (i > 1 || j > 1) {
    if (i > 1 && j > 1 && frechet[i,j] == frechet[i-1,j-1] + distances[i,j]){
      seq1 <- c(i - 1, seq1)
      seq2 <- c(j - 1, seq2)
      i <- i - 1
      j <- j - 1
    } else if (i > 1 && frechet[i,j] == frechet[i-1,j] + distances[i,j]){
      seq1 <- c(i - 1, seq1)
      seq2 <- c(j, seq2)
      i <- i -1
    } else {
      seq1 <- c(i, seq1)
      seq2 <- c(j - 1,seq2)
      j <- j - 1
    }
  }
  
  # returns the coupling, i.e. the ids of points of the first path and the second path
  # the sum of distances between the points is the frechet distance
  return(list("firstseq" = seq1, "secondseq" = seq2))
}

calc_distance <- function(reference_data, test_data, distance_matrix, adjustment){
  # test if tip labels are equal
  if(identical(sort(test_data[[1]]$tip.label), sort(reference_data[[1]]$tip.label))){
    tips <- reference_data[[1]]$tip.label
  } else {
    stop('The tips of the reference and test tree are not equal. The trees cannot be compared. Please check the trees and the tip labels.')
  }
  
  # initialize matrix to save costs for each node
  cost_tree_test <- matrix(data = 0, nrow = (length(test_data[[1]]$tip.label) + length(test_data[[1]]$node.label)), ncol = 2)
  cost_tree_ref <- matrix(data = 0, nrow = (length(reference_data[[1]]$tip.label) + length(reference_data[[1]]$node.label)), ncol = 2)
  
  # loop over all tips
  for (i in 1:length(tips)){
    tip <- tips[i];
    
    # get paths from root to tips
    ref_path <- get_path(reference_data[[1]], tip, reference_data[[2]]);
    test_path <- get_path(test_data[[1]], tip, test_data[[2]]);
    
    # calculate the coupling between the two paths
    # returns a list with paired elements $firstseq and $secondseq including the indices for path1 and path2
    coupling <- frechetCoupling(ref_path$locations, test_path$locations, distance_matrix, sum)
    
    # costs includes the costs for each point in the coupling
    # the node ids for the coupling are given by ref_path$ids[coupling$firstseq] and test_path$ids[coupling$secondseq]
    # assumes symmetrical distances!
    costs <- unname(diag(distance_matrix[ref_path$locations[coupling$firstseq], test_path$locations[coupling$secondseq]]))
    
    # summarize costs per node on the path
    costs_ref <- aggregate(costs ~ ids, data = t(rbind(costs, ids = ref_path$ids[coupling$firstseq])), FUN = sum)
    costs_test <- aggregate( costs ~ ids, data = t(rbind(costs, ids = test_path$ids[coupling$secondseq])), FUN = sum)
    
    # collect costs in matrix; count how often node was in path
    for (p in 1:nrow(costs_test)){
      cost_tree_test[costs_test[p, 1], 1] <- cost_tree_test[costs_test[p, 1], 1] + costs_test[p, 2]
      cost_tree_test[costs_test[p, 1], 2] <- cost_tree_test[costs_test[p, 1], 2] + 1
    }
    
    for (p in 1:nrow(costs_ref)){
      cost_tree_ref[costs_ref[p, 1], 1] <- cost_tree_ref[costs_ref[p, 1], 1] + costs_ref[p, 2]
      cost_tree_ref[costs_ref[p, 1], 2] <- cost_tree_ref[costs_ref[p, 1], 2] + 1
    }
  }
  
  # frechet distance for tree: divide costs per node by number of paths and take sum
  if (adjustment == TRUE){
    dist_test <- sum(cost_tree_test[,1] / cost_tree_test[,2], na.rm = TRUE)
    dist_ref <- sum(cost_tree_ref[,1] / cost_tree_ref[,2], na.rm = TRUE)
  } else {
    dist_test <- sum(cost_tree_test[,1], na.rm = TRUE)
    dist_ref <- sum(cost_tree_ref[,1], na.rm = TRUE)
  }
  # average between distances for ref and test
  return(mean(c(dist_test,dist_ref)))
}

check_dist_matrix <- function(distance_matrix, location_list){
  locations_matrix <- row.names(distance_matrix)
  output <- 1
  # remove duplicates from location_list
  location_list <- unique(location_list)
  # check if each element is in the distance matrix
  for (loc in location_list){
    if (!loc %in% locations_matrix){
      warning(paste("Location", loc, "is not in the distance matrix."), call.=FALSE)
      output <- 0
    }
  }
  return(output)
}

### perform analysis

### load command line arguments ###

library("ape")      # package for phylogenetic trees

args = commandArgs(trailingOnly=TRUE)
reference <- args[1]
test <- args[2]
distances <- args[3]
adjustment <- args[4]

if (is.na(adjustment)){
  adjustment = TRUE
}
if (adjustment != TRUE && adjustment != FALSE){
  stop("Invalid argument. Adjustment must be set to TRUE or FALSE.")
}

### read files ###

reference_data <- read_tree(reference)
test_data <- read_tree(test)

distance_matrix <- as.matrix(read.csv(distances, header = TRUE, row.names = 1, check.names=FALSE, na.strings = ""))
colnames(distance_matrix)[which(is.na(colnames(distance_matrix)))] <- "NA"

### test input ###

# check if all locations are present in the distance matrix
test_loc <- check_dist_matrix(distance_matrix, c(reference_data[[2]]$location, test_data[[2]]$location))

if (test_loc == 0){
  stop("One or more locations are not in the distance matrix. Check your input.")
}

# check if distance matrix is symmetric
symmetric <- isSymmetric(distance_matrix)

if (symmetric == FALSE){
  stop("Distance matrix is not symmetric. Check your input.")
}

# check if sequence names in tree and annotation are equal
if (!identical(sort(c(reference_data[[1]]$tip.label, reference_data[[1]]$node.label)),sort(reference_data[[2]]$label))){
  stop(paste("Sequence ids in tree and annotation are different in files", reference,". Check your input."))
}

if (!identical(sort(c(test_data[[1]]$tip.label, test_data[[1]]$node.label)),sort(test_data[[2]]$label))){
  stop(paste("Sequence ids in tree and annotation are different in files", test,". Check your input."))
}

### calculate distances ###

dist <- calc_distance(reference_data, test_data, distance_matrix, adjustment)
cat(dist) # write output to command line
