library(phangorn)
library(magrittr)

args = commandArgs(trailingOnly=TRUE)
# Read in nexus neighbour net from SplitsTree
Nnet <- read.nexus.networx(args[1])
# Common name prefix for RAxML outputs
tree_prefix <- args[2]
# Path to directory containing input files
inpath <- args[3]
# Output file names
outplot1 = paste(inpath,"/",tree_prefix,"_net.svg",sep="")
outplot2 = paste(inpath,"/",tree_prefix,"_annot_net.svg",sep="")

# Read in RAxML bipartitions
bipartitions <- paste(inpath,"/RAxML_bipartitions.", tree_prefix, sep="")
bootstraps <- paste(inpath,"/RAxML_bootstrap.", tree_prefix, sep="")
raxml.tree <- read.tree(bipartitions)
# Read in RAxML bootstrap trees 
raxml.bootstrap <- read.tree(bootstraps)
# create a vector of labels for the network corresponding to edges in the tree
edge.lab <- createLabel(Nnet, raxml.tree, raxml.tree$edge[,2], "edge")

# find edges that are in the network but not in the tree
edge.col <- createLabel(Nnet, raxml.tree, "black", nomatch="red")

# Confidences values between total, percentage or ratios
x <- addConfidences(Nnet,raxml.tree, scaler = .01)
# Find splits that are in the network but not in the tree
split.col <- createLabel(x, raxml.tree, label="black", "split", nomatch="red")
tipcol <- rep('black', length(x$tip.label))
# Make a vector with our list of populations
# Populations are hard coded for a specific analysis
popIdentifiers <- c("_BS_",
             "_DUN_",
             "_DUR_",
             "_DWESA_",
             "_EBEACH_",
             "_ELDN_",
             "_FRASER_", 
             "_HARVEY_",
             "_HLU_",
             "_HOLE_",
             "_ILUKA_", 
             "_KWIN12",
             "_KWIN18",
             "_LAPER_",
             "_MINNIE_", 
             "_MZN_",
             "_NEW_", 
             "_FAIRE_PISI",
             "_KEMBLA_", 
             "_STJOHN_", 
             "_QMTH_",
             "_STLUC_", 
             "_TMTH_", 
             "_WOL_") 

# Make vector of colors (same order as populations above):
colorsList <- c("sienna", #Boneseed
  "darkgreen", #Dunbogan
  "cornflowerblue", #Durban
  "hotpink", #Dwesa
  "orange", #East Beach
  "lightgoldenrod3", #East London
  "forestgreen", #Fraser Island
  "darkolivegreen", #Harvey Bay
  "lightpink", #Hluleka
  "plum", #Hole in the wall
  "chartreuse4", #Iluka
  "orangered", # Kwinanan 2012
  "orangered1", #Kwinana 2018
  "palegreen3", #La Perouse
  "seagreen", #Minnie Water
  "deepskyblue", #Mtunzini
  "mediumseagreen", #Newcastle
  "sienna3", #Pisifera
  "olivedrab", #Port Kembla
  "palevioletred",# Port St Johns
  "khaki", # Qholora Mouth
  "dodgerblue", #St Lucia
  "lightskyblue", #Tungela Mouth
  "springgreen4" #Wollongong
)

# Set colours for populations
for(i in 1:length(popIdentifiers)){
  tipcol[grep(popIdentifiers[i], x$tip.label)] <- colorsList[i]
}

# Plot without coloured edges/splits
par(mar=c(1,1,1,1))
svg(outplot1, width = 10)
plot(x,"2D",
    tip.color = tipcol, 
    edge.width = .3, 
    cex = .1,
    show.edge.label = FALSE)
dev.off()
# Annotated plot with coloured conflict edges/splits
svg(outplot2, width = 10)
plot(x,"2D",
    tip.color = tipcol, 
    split.color=split.col, 
    edge.width = .3, 
    col.edge.label = "blue", 
    cex = .1,
    show.edge.label = TRUE,
    srt=180)
dev.off()

