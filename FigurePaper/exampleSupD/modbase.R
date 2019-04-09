#!/usr/bin/Rscript
## EXAMPLE D
## Description: Get lists of PDBs containing a 3 RNA modified residues,
##              then filter by resolution and execute DSSR. Plot motif in which
##              the modified residues are found.  
## Author: Diego Gallego

## ----------------------------------------------------------------------------
## Preprocess
## ----------------------------------------------------------------------------
## Load necessary packages for the example
library(veriNA3d)
library(ggplot2) ## Specific dependency for this example!
library(gridExtra) ## Specific dependency for this example!
library(scales) ## Specific dependency for this example!

## Load specific functions for the example
source("dependencies.R")

## ----------------------------------------------------------------------------
## Get Raw data
cat("Getting and cleaning data\n")
## ----------------------------------------------------------------------------
modres <- c("1MA", "6MA", "4SU")
motifs <- c("hairpin-loop", "junction-loop", "ss-non-loop", "internal-loop")
motifs2 <- c("helix")#, "helix-end")

out <- lapply(modres, function(res) {
    ## Query which structures contain the modres of interest
    pdblist <- toupper(unlist(queryAPI(ID=res, API="ebi", verbose=TRUE,
                                string1="pdb/compound/in_pdb/", string2="")))

    ## Query resolution and keep only structures under 3 A
    resol <- applyToPDB(pdblist, queryResol)
    pdblist <- resol[resol[, 2] < 3, 1]

    output <- c(); counter <- 0
    ## Loop over list of PDB IDs
    for (i in seq_along(pdblist)) {
        tryCatch({
                ## Execute DSSR
                dssrdata <- dssr(pdblist[i])
## ----------------------------------------------------------------------------
## Clean data. Raw -> Tidy
## ----------------------------------------------------------------------------
                ## Extract data of interest
                data <- dssrdata$models$parameters$nts[[1]]
                summary <- data$summary[data$nt_name == res]
                for (j in seq_along(summary)) {
                    dat <- unlist(strsplit(summary[j], ","))
                    ind <- which(dat %in% motifs)
                    if (length(ind) == 0) {
                        ind <- which(dat %in% motifs2)
                    }
                    if (length(ind) != 0) {
                        counter <- counter + 1
                        output[counter] <- dat[ind]
                        in2 <- which(data$nt_name == res)
                        ## Save neigbouring sequence
                        write(c(dat[ind], 
                                paste(data$nt_name[(in2-2):(in2+2)], collapse="-")), 
                                file="seq.txt", ncolumns=2, append=T, sep=" ")
                    }
                }
            ## If DSSR fails with the given structure, do nothing
            }, error=function(e) {})
    }
    return(output)
})
## ----------------------------------------------------------------------------
## Exploratory analysis. Plot
cat("Plotting data\n")
## ----------------------------------------------------------------------------
p <- lapply(1:3, function(i) {
    plotM(out[[i]], res=modres[i], motifs=motifs, motifs2=motifs2)
})
tiff("exampleD.tiff", height=10.5, width=24.75, units="cm", res=300)
grid.arrange(p[[1]], p[[2]], p[[3]], ncol=3)
dev.off()
