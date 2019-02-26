## Load data and library
load("PO.RData")
require(pROC)

## AUC calculation
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

## AUC for stanford11 
AUCs_stanford <- rep(0,12,1)
AUCs_stanford_ci <- matrix(0,12,3)

tm <- AUCcalc(GSE28750_expr, GSE28750_sample, stanford11.genes)
AUCs_stanford[1] <- tm$auc
AUCs_stanford_ci[1,] <- t(as.matrix(tm$ci))
tm <- AUCcalc(GSE32707_expr, GSE32707_sample, stanford11.genes)
AUCs_stanford[2] <- tm$auc
AUCs_stanford_ci[2,] <- t(as.matrix(tm$ci))
tm <- AUCcalc(GSE40012_expr, GSE40012_sample, stanford11.genes)
AUCs_stanford[3] <- tm$auc
AUCs_stanford_ci[3,] <- t(as.matrix(tm$ci))
tm <- AUCcalc(GSE66099_expr, GSE66099_sample, stanford11.genes)
AUCs_stanford[4] <- tm$auc
AUCs_stanford_ci[4,] <- t(as.matrix(tm$ci))
tm <- AUCcalc(GSE36809_expr1, GSE36809_sample1, stanford11.genes)
AUCs_stanford[5] <- tm$auc
AUCs_stanford_ci[5,] <- t(as.matrix(tm$ci))
tm <- AUCcalc(GSE36809_expr2, GSE36809_sample2, stanford11.genes)
AUCs_stanford[6] <- tm$auc
AUCs_stanford_ci[6,] <- t(as.matrix(tm$ci))
tm <- AUCcalc(GSE36809_expr3, GSE36809_sample3, stanford11.genes)
AUCs_stanford[7] <- tm$auc
AUCs_stanford_ci[7,] <- t(as.matrix(tm$ci))
tm <- AUCcalc(GSE36809_expr4, GSE36809_sample4, stanford11.genes)
AUCs_stanford[8] <- tm$auc
AUCs_stanford_ci[8,] <- t(as.matrix(tm$ci))
tm <- AUCcalc(GSE36809_expr5, GSE36809_sample5, stanford11.genes)
AUCs_stanford[9] <- tm$auc
AUCs_stanford_ci[9,] <- t(as.matrix(tm$ci))
tm <- AUCcalc(GSE65682_expr, GSE65682_sample, stanford11.genes)
AUCs_stanford[10] <- tm$auc
AUCs_stanford_ci[10,] <- t(as.matrix(tm$ci))
tm <- AUCcalc(GSE74224_expr, GSE74224_sample, stanford11.genes)
AUCs_stanford[11] <- tm$auc
AUCs_stanford_ci[11,] <- t(as.matrix(tm$ci))
tm <- AUCcalc(EMEXP3589_expr, EMEXP3589_sample, stanford11.genes)
AUCs_stanford[12] <- tm$auc
AUCs_stanford_ci[12,] <- t(as.matrix(tm$ci))
colnames(AUCs_stanford_ci) <- c("95% lower", "AUC", "95% upper")
