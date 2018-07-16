#!/usr/bin/env Rscript

library("ggplot2")
library("MASS")

# visualize pairwise Frechét distances using multidimensional scaling

args = commandArgs(trailingOnly=TRUE)
distancefile <- args[1]

values <- read.table(distancefile, stringsAsFactors=FALSE)

objects <- unique(c(values$V1, values$V2))

# initialize matrix
values.matrix <- matrix(data=0, nrow=length(objects), ncol=length(objects))
rownames(values.matrix) <- objects
colnames(values.matrix) <- objects

# fill matrix
for (i in 1:nrow(values)){
  values.matrix[values[i,1], values[i,2]] <- values[i,3]
  values.matrix[values[i,2], values[i,1]] <- values[i,3]
}

# perform metric multidimensional scaling (principal component analysis)
fit_metric <- cmdscale(values.matrix, eig = TRUE, k = 2)
mds_data <- data.frame(x = fit_metric$points[, 1], y = fit_metric$points[, 2], label = rownames(fit_metric$points))

# nonmetric multidimensional scaling
fit_nonmetric <- isoMDS(values.matrix, y = cmdscale(values.matrix, k = 2), k = 2)
mds_data <- data.frame(x = fit_nonmetric$points[, 1], y = fit_nonmetric$points[, 2], label = rownames(fit_nonmetric$points))

png(file="pairwiseFrechetDistances.png", width=8, height=3, units="in", res=300)

ggplot(mds_data, aes(x,y)) + 
  geom_point() +
  geom_text(aes(x=x-1000, y=y+500, label=label)) +
  coord_fixed(ratio=1)

#ggsave("Distances.png", width = 8, height = 3)
