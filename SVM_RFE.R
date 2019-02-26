#########################################################################
##
##  SVM_RFE.R: Performing Support Vector Machine with Recursive 
##  Feature Elimination (SVM-RFE) for the co-normalized microarray data 
##  
#########################################################################


## 1. Start from the results generated from COCONUT.R

## 2. Generation of normalized expression matrix for Stanford82
tid <- which(is.element(row.names(GSEs.COCO.combined$genes), stanford82.genes) == 1)
normData <- GSEs.COCO.combined$genes[tid, ]
sampleinfo <- GSEs.COCO.combined$pheno

## 3. feature selection
## run pathClass
library(pathClass)
library(mdatools)

x <- t(normData)
xAuto <- prep.autoscale(x, center = T, scale = T)
y <- as.factor(sampleinfo$cond)

res.rfe <- crossval(xAuto, y, DEBUG=TRUE, theta.fit=fit.rfe, folds=5, repeats=100, parallel=TRUE, Cs=10^(-3:3))
extractFeatures(res.rfe, toFile=TRUE, fName="SVMRFE_frequency_Auto_new_100.csv")

## AUC calculation by increasing the number of genes sorted by frequency
data <- read.table("SVMRFE_frequency_Auto_100.csv", header=T, sep=",")
geneOrder <- data$time.choosen.Var1

set.seed(7)
library(caret)
library(e1071)
library(ROCR)

X <- as.data.frame(t(normData))
genes <- colnames(X)
y <- sampleinfo$cond
y[which(sampleinfo$cond == 1)] = 'Sepsis'
y[which(sampleinfo$cond == 0)] = 'SIRS'

nRow <- nrow(X)
nTimes <- 100
kFolds <- 5
n <- 82
AUCvals <- matrix(0, n, 1)

for ( k in 1:n ) {
  print(sprintf('%d', k))
  if (k == 1) {
    features <- geneOrder[1]
    sids <- which(is.element(genes, features) == 1)
    sGenes <- genes[sids]
    Xsub <- as.vector(X[, sids])
    AUCvals[k] = auc(roc(y, Xsub, ci=T))
  } else {
    features <- geneOrder[1:k]
    sids <- which(is.element(genes, features) == 1)
    sGenes <- genes[sids]
    Xsub <- as.data.frame(prep.autoscale(X[, sids], center=T, scale = T))
    Xsub <- cbind(Xsub, y)
    colnames(Xsub) <- c(sGenes, 'y')
    prob <- matrix(0, nRow, 1)
    AUCtm <- matrix(0, nTimes, 1)
    for ( i in 1:nTimes ) {
      set.seed(i)
      fids <- createFolds(y, k = kFolds, list = T, returnTrain = F)
      for ( j in 1:kFolds ) {
        testid <- fids[j]$Fold
        testing <- as.data.frame(Xsub[testid, ])
        training <- as.data.frame(Xsub[-testid, ])
        svmModel <- svm(y ~ ., data = training, kernel = 'linear', probability = T)
        pred <- predict(svmModel, testing, probability = T)
        prob[testid] <- as.matrix(attr(pred, 'probabilities')[, 2])
      }
      AUCtm[i] <- auc(roc(y, as.vector(prob), ci=T))
    }
    AUCvals[k] <- mean(AUCtm)
  }
}

row.names(AUCvals) <- geneOrder
maxid <- order(AUCvals, decreasing=T)[1]
RFEgenes <- row.names(AUCvals)[1:maxid]