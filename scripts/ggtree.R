#!/usr/bin/env Rscript

library("ape")
library("Biostrings")
library("ggplot2")
library("ggtree")
library("colorspace")

args = commandArgs(trailingOnly=TRUE)
intree <- args[1]
popmap <- args[2]

outname1 <- paste(intree,"_nobrln_tree.svg",sep="") 
outname2 <- paste(intree,"_brln_tree.svg",sep="") 

pops <- read.csv(popmap,header = F, sep = "\t")
# Get number of populations
numpop <- length(unique(pops$V2))

#Define OTUs
poplist <- split(pops$V1, pops$V2)


raxml <- read.tree(intree)
#Assign OTU
mytree <- groupOTU(raxml, poplist)
svg(outname1, width = 10)
ggtree(mytree, aes(color=group),branch.length="none", layout="rectangular") + geom_tiplab(size=1) +
  scale_color_manual(values=c("black", rainbow_hcl(numpop))) + theme(legend.position="right")+ 
  geom_text2(size = 2,aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 70))+
  scale_size_manual(values=c(1, .1))  + ggplot2::xlim(0, 60)
dev.off()
svg(outname2, width = 10)
ggtree(mytree, aes(color=group), layout="rectangular") + geom_tiplab(size=1) +
  scale_color_manual(values=c("black", rainbow_hcl(numpop)))  + theme(legend.position="right")+ 
  geom_text2(size = 2,aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 70))+
  scale_size_manual(values=c(1, .1)) 
dev.off()
