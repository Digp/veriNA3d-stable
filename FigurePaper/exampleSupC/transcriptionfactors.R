#!/usr/bin/Rscript
## EXAMPLE C
## Description: From a dataset of DNA-TF, plot the most frequent amino acids
##              interacting with phosphate, sugar and base.
## Author: Diego Gallego

## ----------------------------------------------------------------------------
## Preprocess
## ----------------------------------------------------------------------------
## Load necessary packages for the example
library(veriNA3d)
library(bio3d)

plotAA <- function(ta, main) {
    color <- rainbow(n=20)
    names(color) <- c("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", 
                        "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
                        "THR", "TRP", "TYR", "VAL")
    par(mar=c(0, 3, 1, 0.6))
    a <- barplot(ta[5:1]/sum(ta), xlim=c(0, 0.7), las=2, border=0,
                    #cex.lab=1, cex.axis=1, font=1, col=color[names(ta[1:5])],
                    xaxt='n', font=1, col=color[names(ta[5:1])],
                    cex.names=1.5, main=main, cex.main=1.5, horiz=TRUE)
}

## ----------------------------------------------------------------------------
## Get Raw data
cat("Getting data\n")
## ----------------------------------------------------------------------------

path <- "transcription_factor_pdb/"
pdblist <- dir(path, pattern=".ent.gz")
pdblist <- gsub(".ent.gz", "", pdblist)

## Get structural data
aantinfo <- pipeProtNucData(pdbID=pdblist, path=path, 
                            extension=".ent.gz", cutoff=10, cores=1)

## ----------------------------------------------------------------------------
## Clean data. Raw -> Tidy
cat("Cleaning data\n")
## ----------------------------------------------------------------------------

## Filter to keep only closest contacts using distance < 3 Angstroms
intdata <- aantinfo[aantinfo$distance < 3, ]

## Create objects of interest
P <- c("OP1", "OP2", "OP3")
S <- c("O4'", "O3'", "O5'")
B <- c("N1", "N2", "N3", "N4", "N6", "N7", "O2", "O4", "O6")

ta_P <- sort(table(intdata[intdata$elety_A %in% P, "resid_B"]), decreasing=TRUE)
ta_S <- sort(table(intdata[intdata$elety_A %in% S, "resid_B"]), decreasing=TRUE)
ta_B <- sort(table(intdata[intdata$elety_A %in% B, "resid_B"]), decreasing=TRUE)

## ----------------------------------------------------------------------------
## Exploratory analysis. Plot
cat("Plotting data\n")
## ----------------------------------------------------------------------------

total <- sum(ta_P, ta_S, ta_B)
tiff("exampleC.tiff", height=10.5, width=16.5, units="cm", res=300)
par(mfrow=c(3, 1), oma=c(4.5, 1, 1, 1))
## Phosphate interactions plot
plotAA(ta_P, main=paste("Phosphate ", round(sum(ta_P)/total*100), "%", sep=""))
## Sugar interactions plot
plotAA(ta_S, main=paste("Sugar ", round(sum(ta_S)/total*100), "%", sep=""))
## Base interactions plot
plotAA(ta_B, main=paste("Base ", round(sum(ta_B)/total*100), "%", sep=""))
axis(side=1, lwd=1.5, cex.axis=2)
title(xlab="Density", cex.lab=2, outer=TRUE)
dev.off()
