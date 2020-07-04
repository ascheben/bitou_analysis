#!/usr/bin/env Rscript

library(pophelper)
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)
inpath <- args[1]
labels <- args[2]
grplabels <- args[3]
lenK <- args[4]

shiny <- c("#1D72F5","#DF0101","#77CE61", "#FF9326","#A945FF","#0089B2","#FDF060","#FFA6B2","#BFF217","#60D5FD","#CC1577","#F2B950","#7FB21D","#EC496F","#326397","#B26314","#027368","#A4A4A4","#610B5E")


#Add all files in directory to list
sfiles <- list.files(path = inpath, pattern = "\\.meanQ$", full.names=TRUE)
#read files in from list and add individual labes from file
#files can be modifief to show cultivar names
slist <- readQ(files=sfiles,indlabfromfile=T)
#read individual labels
inds <- read.delim(labels,header=F,stringsAsFactors=FALSE)
#add ind names as rownames to all tables
if(length(unique(sapply(slist,nrow)))==1) slist <- lapply(slist,"rownames<-",inds$V1)


grouplabset <- read.delim(grplabels, header=T,stringsAsFactors=F)


out_all <- paste(inpath,"/structure_plot_lenK",lenK,".pdf",sep="")
#plot selected K with inidividuals in groups
allk <- plotQ(slist[1:lenK],imgoutput="join",returnplot=T,exportplot=F, basesize=11,
              quiet=T, showticks = T,grplab=grouplabset,ordergrp=TRUE,showindlab=T,useindlab = T,
              indlabsize = 4,grplabsize=4,splabsize=5, clustercol=shiny)
pdf(out_all)
grid.arrange(allk$plot[[1]])
dev.off()

#loop to make plot for each K
for (i in 1:lenK) {
        single_out <- paste(inpath,"/structure_plot_K",i,".pdf",sep="");
        singlek <- plotQ(slist[i],returnplot=T,exportplot=F, basesize=11,
            quiet=T, showticks = T,grplab=grouplabset,ordergrp=TRUE,showindlab=T,useindlab = T,
            indlabsize = 4,grplabsize=4,splabsize=5, clustercol=shiny);
        pdf(single_out, width = 12, height = 5);    
        grid.arrange(singlek$plot[[1]])
        dev.off();
        }
