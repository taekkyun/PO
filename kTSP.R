#########################################################################
##
##  kTSP.R: Performing k-Top Scoring Pairs classifier (kTSP) using 
##  switchbox R package
##  
#########################################################################

## Load data and library
load("PO.RData")
require(switchBox)

## k-Top scoring pair(s)
matTraining <- cbind(GSE28750_expr,
                    GSE32707_expr,
                    GSE40012_expr,
                    GSE66099_expr,
                    GSE36809_expr1,
                    GSE36809_expr2,
                    GSE36809_expr3,
                    GSE36809_expr4,
                    GSE36809_expr5)
matTraining_group <- as.factor(c(as.character(GSE28750_sample$group),
                         as.character(GSE32707_sample$group),
                         as.character(GSE40012_sample$group),
                         as.character(GSE66099_sample$group),
                         as.character(GSE36809_sample1$group),
                         as.character(GSE36809_sample2$group),
                         as.character(GSE36809_sample3$group),
                         as.character(GSE36809_sample4$group),
                         as.character(GSE36809_sample5$group)))
matTraining_group_number <- as.numeric(matTraining_group) - 1
row.names(matTraining) <- stanford82.genes

classifier <- SWAP.KTSP.Train(matTraining, matTraining_group, krange=c(3:15))
ktspStatDefault <- SWAP.KTSP.Statistics(inputMat = matTraining, classifier = classifier)
testPrediction <- SWAP.KTSP.Classify(matTraining, classifier)

## AUC calculation
## training using optimal set of TSPs
ktsp <- KTSP.Train(matTraining, matTraining_group_number, n=3) # n was decided during previous step
AUCsTSP <- rep(0, 12)
pVal <- matrix(0, 12, 1)

## testing for GSE28750
i <- 1
matTesting <- GSE28750_expr
row.names(matTesting) <- stanford82.genes
matTesting_group_number <- as.numeric(factor(GSE28750_sample$group)) - 1

pred <- KTSP.Classify(matTesting, ktsp, combineFunc = sum)
perf <- roc(!matTesting_group_number, pred, ci=FALSE, plot=FALSE)
AUCsTSP[i] <- auc(perf)

tROC1 <- AUCcalc(GSE28750_expr, GSE28750_sample, stanford11.genes)
tROC2 <- roc(!matTesting_group_number, pred, ci=FALSE, plot=FALSE)
ROCres <- roc.test(tROC1, tROC2, method="delong")
pVal[i,1] <- ROCres$p.value

## testing for Other datasets using the same code