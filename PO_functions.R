# AUC calculation
AUCcalc <- function (data, sample, genes) {
  data_expr <- data
  data_sample <- sample
  
  sepsisid <- which(data_sample$group == "Sepsis")
  sirsid <- which(data_sample$group == "SIRS")
  gid <- c()
  for ( i in 1:length(genes) ){
    tid <- which(stanford82.genes == genes[i])
    gid <- c(gid, tid)
  }
  folddirection <- stanford82.foldDirection[gid]
  
  tMat <- data_expr[gid, c(sepsisid,sirsid)]
  
  # check zero in matrix
  tid <- which(tMat[,1] != 0)
  tMat <- tMat[tid, ]
  folddirection <- folddirection[tid]
  
  # AUC calculation
  numcol <- ncol(tMat)
  posid <- which(folddirection > 0)
  negid <- which(folddirection < 0)
  
  SMS <- rep(0,numcol)
  for ( i in 1:numcol ) {
    if ( length(posid) > 0 && length(negid) > 0 ) {
      posExp <- exp(mean(log(tMat[posid,i])))
      negExp <- exp(mean(log(tMat[negid,i])))
      weight <- length(negid) / length(posid)
      SMS[i] <- posExp - (negExp * weight)
    } else if ( length(posid) > 0 && length(negid) == 0 ) {
      posExp <- exp(mean(log(tMat[posid,i])))
      SMS[i] <- posExp
    } else if ( length(posid) == 0 && length(negid) > 0 ) {
      negExp <- exp(mean(log(tMat[negid,i])))
      SMS[i] <- -negExp
    } else {
      SMS[i] <- 0
    }
  }
  
  class <- factor(data_sample$group[c(sepsisid,sirsid)], levels = c("Sepsis", "SIRS"))
  library(pROC)
  roc_obj <- roc(class, SMS, ci=T)
  return(roc_obj)
  
}

# P-value calculation
ROCcomparePval <- function (data, sample, AUCtable) {
  sid <- c(1, 2, 3, 4, 5, 9)
  pvalueMat <- rep(1, 6)
  for ( i in 1:length(sid) ) {
    tAUCs <- cbind(substitutingGenesMatrix[,i], substitutedGenesMatrix[,i], AUCtable[,i])
    tAUCs <- tAUCs[which(is.na(tAUCs[,1]) == 0), ]
    tid <- which(tAUCs[, 3] == median(as.numeric(tAUCs[, 3])))
    ttGenes <- tAUCs[tid[1], 1:2]
    j <- which(stanford11.genes == ttGenes[1])
    tGenes <- stanford11.genes
    tGenes[j] <- ttGenes[2]
    tm <- AUCcalc(data, sample, stanford11.genes)
    tROC1 <- tm$auc
    tm <- AUCcalc(data, sample, stanford11.genes)
    tROC2 <- tm$auc
    ROCres <- roc.test(tROC1, tROC2, method="delong")
    pvalueMat[i] <- ROCres$p.value
  }
  return(pvalueMat)
  
}

# Comparison of ROCs
ROCcompare6GenesPval <- function (data, sample, AUCtable) {
  tAUCs <- cbind(tGeneList[selid,], AUCtable[selid])
  tid <- which(tAUCs[, 12] == median(as.numeric(tAUCs[, 12])))
  tGenes <- tAUCs[tid[1], 1:11]
  tROC1 <- AUCcalc(data, sample, stanford11.genes)
  tROC2 <- AUCcalc(data, sample, tGenes)
  ROCres <- roc.test(tROC1, tROC2, method="delong")
  return(ROCres$p.value)
}
