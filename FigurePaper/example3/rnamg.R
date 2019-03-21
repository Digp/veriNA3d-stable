#!/usr/bin/Rscript
## EXAMPLE 3
## Description: From the whole PDB, create a dataset of RNA structures with Mg
##              and a threshold resolution of 2 Angstroms. Then, measure and
##              plot the distances between Mg and other atom types.
## Author: Diego Gallego

## ----------------------------------------------------------------------------
## Preprocess
## ----------------------------------------------------------------------------
## Load necessary packages for the example
library(veriNA3d)
library(ggplot2) ## Specific dependency for this example!
library(reshape2) ## Specific dependency for this example!

## ----------------------------------------------------------------------------
## Get Raw data
cat("Getting data\n")
## ----------------------------------------------------------------------------
## From the whole PDB, prepare dataset of interest ============================
## Retrieve ALL entries in the PDB (including latest entries)
pdblist <- queryEntryList(justIDs=FALSE)

## Keep only Nucleic Acid structures solved by diffraction techniques
nalist <- pdblist[which(grepl("nuc", pdblist$type, perl=TRUE) & 
                        pdblist$technique == "diffraction"), "pdbID"]

## Query which structures contain MG
mg <- toupper(unlist(queryAPI(ID="MG", API="ebi", verbose=TRUE,
                                string1="pdb/compound/in_pdb/", string2="")))
mgnalist <- nalist[which(nalist %in% mg)]

## Query resolution and keep only structures under 2 A
resol <- applyToPDB(mgnalist, queryResol)
mgnalist <- resol[resol[, 2] < 2, 1]

## Check entities in the structures and keep only PDBs with RNA
naent <- applyToPDB(mgnalist, countEntities, as.df=FALSE)
naent <- as.data.frame(do.call(rbind, naent), stringsAsFactors=FALSE)
mgrnalist <- mgnalist[which(naent$RNA > 0)]
## ============================================================================

## Loop over list of PDB IDs
distances <- do.call(rbind, lapply(seq_along(mgrnalist), function(i) {
    ## Download mmCIF and get entity of Mg
    cif <- cifParser(mgrnalist[i])
    MG <- cifEntity(cif)[grepl("MAGNE", cifEntity(cif)$pdbx_description), "id"]

    ## Measure distances between Mg and other entities
    return(measureEntityDist(cif, refent=MG, entities="all", cutoff=4, n=NULL))
}))

## ----------------------------------------------------------------------------
## Clean data. Raw -> Tidy
cat("Cleaning data\n")
## ----------------------------------------------------------------------------
Ow <- distances$distance[distances$resid_B == "HOH" & distances$elety_B == "O"]
distances <- distances[distances$resid_B %in% c("A", "C", "G", "U"), ]
Ob <- distances$distance[distances$elety_B %in% c("O2", "O4", "O6")]
Or <- distances$distance[distances$elety_B %in% c("O2'", "O3'", "O4'", "O5'")]
Op <- distances$distance[distances$elety_B %in% c("OP1", "OP2", "OP3")]
Nb <- distances$distance[distances$elety_B %in% c("N1", "N2", "N3", "N7")]
## Data frame with distances per type of aotm
dist <- list(Nb=Nb, Ob=Ob, Or=Or, Op=Op, Ow=Ow)
df <- melt(dist, value.name="Distance"); names(df)[2] <- "atom"
df$atom <- factor(df$atom, levels=c("Ow", "Op", "Or", "Ob", "Nb"))

## ----------------------------------------------------------------------------
## Exploratory analysis. Plot
cat("Plotting data\n")
## ----------------------------------------------------------------------------
tiff("example3.tiff", height=10.5, width=16.5, units="cm", res=300)
colors <- c("#87CEEB", "#6A5ACD", "#F08080", "#DC143C", "#4169E1")
ggplot(df, aes(x=Distance, fill=atom)) + ylab("Frequency") + 
        scale_fill_manual(values=colors) + scale_y_continuous(breaks=NULL) +
        scale_x_continuous(breaks=seq(1.2, 4, 0.4)) +
        geom_histogram(binwidth=0.03, color="white") +
        theme(panel.grid.minor=element_blank())
dev.off()
