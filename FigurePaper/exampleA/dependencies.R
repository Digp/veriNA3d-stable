## "cache" saves trinucleotide backbones in RAM for latter RMSD calculations
cache <-
function(ntID, name, ntinfo) {
    ## Trim trinucleotide backbone
    pdb <- trimByID(ntID=ntID, ntinfo=ntinfo, prev=1, post=1,
                    alt=c("uniq"), cutoff=0, justbackbone=TRUE,
                    verbose=FALSE)
    ## Save as an independent pdb file in RAM
    assign(name, pdb, envir=.GlobalEnv)
}

## "pwrmsd" calculates pair wise RMSD between nucleotides presaved in RAM
pwrmsd <-
function(ntIDs) {
    pdb1 <- get(paste("nt", ntIDs[1], sep=""))
    pdb2 <- get(paste("nt", ntIDs[2], sep=""))
    sel1 <- atom.select(pdb1)
    sel2 <- atom.select(pdb2)
    fit=fit.xyz(pdb1$xyz, pdb2$xyz,
                fixed.inds=sel1$xyz, mobile.inds=sel2$xyz)
    return(rmsd(pdb1$xyz[,sel1$xyz],fit[,sel2$xyz]))
}

## "filter_pyle" reproduces the iterative process of RMSD comparison to chose
## an helical nucleotide of reference
filter_pyle <- 
function(ntinfo, helicalntID, RMSD, bandwidths=c(40,40), cutoff=0.85) {

    twideMatrix <- t(acast(RMSD, nt1~nt2, value.var="rmsd"))
    twideMatrix<-rbind(rep(NA,length(helicalntID)-1),twideMatrix)
    twideMatrix<-cbind(twideMatrix,rep(NA,length(helicalntID)))
    rownames(twideMatrix)[1]<-helicalntID[1]
    colnames(twideMatrix)[length(helicalntID)]<-helicalntID[length(helicalntID)]

    Matrix <- twideMatrix
    while (length(helicalntID) > 1) {
        Matrix <- Matrix[which(rownames(Matrix) %in% helicalntID),
                            which(colnames(Matrix) %in% helicalntID)]

        column_sums <- colSums(Matrix, na.rm=TRUE)
        row_sums <- rowSums(Matrix, na.rm=TRUE)
        RMSDaverages <- (column_sums + row_sums)/(length(row_sums) - 1)
        index <- which.min(RMSDaverages)
        lowestaverage <- RMSDaverages[index]
        reference_nt <- helicalntID[index]
        A <- Matrix[which(rownames(Matrix) == reference_nt), ]
        names(A) <- colnames(Matrix)
        A <- A[complete.cases(A)]
        B <- Matrix[, which(colnames(Matrix) == reference_nt)]
        names(B) <- colnames(Matrix)
        B <- B[complete.cases(B)]
        reference_nt_RMSD <- append(A,B); rm(list=c("A","B"))
        helicalntID <- sort(append(reference_nt, as.numeric(
                            names(reference_nt_RMSD)[which(
                            reference_nt_RMSD < lowestaverage)])))
    }

    return(helicalntID)
}
