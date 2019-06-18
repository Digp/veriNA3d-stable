#!/usr/bin/Rscript
## EXAMPLE A
## Description: From Leontis non-redundant RNA list, get eta-theta plot
##              following Pyle's method.
## Author: Diego Gallego

## ----------------------------------------------------------------------------
## Preprocess
## ----------------------------------------------------------------------------
## Load necessary packages for the example
library(veriNA3d)
library(bio3d)
library(reshape2) ## Specific dependency for this example!

## Load specific functions for the example
source("dependencies.R")

## ----------------------------------------------------------------------------
## Get Raw data
cat("Getting data\n")
## ----------------------------------------------------------------------------

## Get latest Leontist non-redundant list of RNA structures
infochains <- getRNAList(threshold="2.5A", as.df=TRUE)

## Get structural data
ntinfo <- pipeNucData(pdbID=infochains$pdb, model=infochains$model,
                        chain=infochains$chain, cores=1)

## ----------------------------------------------------------------------------
## Clean data. Raw -> Tidy
cat("Cleaning data\n")
## ----------------------------------------------------------------------------

## Subset north nucleotides
north <- ntinfo[cleanByPucker(ntinfo, surenorth=TRUE), ]

## Find non-broken canonical nucleotides
id <- which(north$puc_valid==TRUE & north$Break==FALSE & north$big_b == FALSE & 
            north$eta_valid == TRUE & north$theta_valid == TRUE &
            grepl("^[GAUC]-[GAUC]-[GAUC]$", north$localenv, perl=TRUE))
usefulnt <- north$ntID[id]

## Filter helical nucleotides before ploting =================================

## Create all the necessary trinucleotides for pair-wise RMSD calculations
invisible(mapply(FUN=cache, ntID=usefulnt, name=paste("nt", usefulnt, sep=""),
                    MoreArgs=list(ntinfo=ntinfo)))

## Find helical nucleotides
HDR <- findHDR(ntID=usefulnt, ntinfo=ntinfo, x="eta", y="theta", SD_DENS=15)

## Compute pair-wise RMSD between helical nucleotides
pairwise <- t(combn(HDR[[1]], 2))
RMSD_calc <- apply(FUN=pwrmsd, MARGIN=1, X=pairwise)
RMSD <- data.frame(nt1=pairwise[, 1], nt2=pairwise[, 2], rmsd=RMSD_calc)

## Apply Pyle method
ref <- filter_helical(ntinfo=ntinfo, helicalntID=HDR[[1]], RMSD=RMSD)

pairs <- data.frame(usefulnt=usefulnt, reference=rep(ref, length(usefulnt)))
RMSD2 <- apply(FUN=pwrmsd, MARGIN=1, X=pairs)

## ===========================================================================
## ----------------------------------------------------------------------------
## Exploratory analysis. Plot
cat("Plotting data\n")
## ----------------------------------------------------------------------------
tiff("exampleA.tiff", height=10.5, width=16.5, units="cm", res=300)
plot2D(ntinfo=ntinfo, x="eta", y="theta", etatheta=TRUE,
        ntID=pairs[which(RMSD2 > 0.85), 1])
dev.off()
