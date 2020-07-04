#!/usr/bin/env Rscript
library(phangorn)
library(magrittr)

# AIM
# The aim of this script is to explore a SplitsTree network and output plots

# PREPARATION
# An IBS distance matrix was calculated from a VCF using TASSEL 5.2.50 with default settings
# The comment rows starting with "#" were deleted from the file, and the file suffix set to ".dst"
# This allowed the files to be loaded into SplitsTree 4.14.8 where a NeighborNet was calculated
# The resulting net was saved as a nexus file and loaded into R



## RAxML best-known tree with bipartition support (from previous analysis)
raxml.tree <- read.tree("RAxML_bipartitions.bitou_core_m2_rename")
## RAxML bootstrap trees (from previous analysis)
raxml.bootstrap <- read.tree("RAxML_bootstrap.bitou_core_m2.rename.nohomhet.rename.tre")
Nnet <- read.nexus.networx("bitou_core_m2_ibs.dst.nex")
#par(mfrow=c(1,2), mar=c(1,1,1,1)) # Setting plot parameters
### Plotting trees with support values:
##  RAxML
#plot(raxml.tree, cex = 0.3)
#nodelabels(raxml.tree$node.label, adj = c(1, 0), frame = "none", cex = 0.3)

#par(mfrow=c(1,1)) # Setting plot parameters
# NeighbourNet
#plot(Nnet,"2D",cex = 0.3)

# create a vector of labels for the network corresponding to edges in the tree
edge.lab <- createLabel(Nnet, raxml.tree, raxml.tree$edge[,2], "edge")
# could be also 1:27 instead of raxml.tree$edge[,2]

# Show the correspondingly labelled tree and network in R
#par(mfrow=c(1,2))  
#plot(raxml.tree, "u", rotate.tree = 180, cex=.3) 
#edgelabels(raxml.tree$edge[,2],col="blue", frame="none", cex=.7)

# find edges that are in the network but not in the tree
#edge.col <- rep("black", nrow(Nnet$edge))
#edge.col[ is.na(edge.lab) ] <- "red"
# or a simpler alternative...
edge.col <- createLabel(Nnet, raxml.tree, "black", nomatch="red")

#x <- plot(Nnet, edge.label = edge.lab, show.edge.label = T, "2D", edge.color = edge.col,
#          col.edge.label = "blue", cex=.7)

# the scaler argument multiplies the confidence values. This is useful to switch
# confidences values between total, percentage or ratios.   
x <- addConfidences(Nnet,raxml.tree, scaler = .01)
# find splits that are in the network but not in the tree
#split.col <- rep("black", length(x$splits))
#split.col[ !matchSplits(as.splits(x), as.splits(raxml.tree)) ] <- "red"

# simpler alternative...
split.col <- createLabel(x, raxml.tree, label="black", "split", nomatch="red")
# Plotting in R

tipcol <- rep('black', length(x$tip.label))
# make a vector with populations
Pops <- c("WAU_", "_QLD_", "_NSW_", "RSA_", "_PISI_", "_BS_")
#Pops <- c("WAU_", "_QLD_", "_NSW_", "RSA_", "_PISI_")
# make a vector of color we want:
colorsList <-c("pink", "darkolivegreen3", "green", "orange", "yellow","purple")

# replace colours where grep gives "TK" as red, etc in a loop
for(i in 1:length(Pops)){
  tipcol[grep(Pops[i], x$tip.label)] <- colorsList[i]
}

par(mar=c(1,1,1,1))
out.x <- plot(x,"2D", show.edge.label=TRUE, tip.color = tipcol, split.color=split.col, edge.width = .3, col.edge.label = "blue", cex = .3)
legendcol <- gsub("black", "darkolivegreen3", tipcol)
#Pops <- c("WAU_", "_QLD_", "_NSW_","_PISI_","RSA_")
legendtext <- gsub("_", "", Pops)
legend("topright",legendtext,fill=unique(legendcol), cex = 0.5)

write.nexus.networx(out.x,"RAxML_bootstrap.bitou_core_m2_consensus_network_015.new.nxs")

#y <- addConfidences(Nnet, as.splits(raxml.bootstrap))
#edge.col <- createLabel(y, raxml.tree, label="black", "edge", nomatch="grey")
#y <- plot(y,"2D",show.edge.label=TRUE,srt = 180, edge.color=edge.col,edge.width = .3, cex = .2)
