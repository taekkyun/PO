#########################################################################
##
##  COCONUT.R: Co-normalization procedure for pooled analysis 
##  of microarray data by batch-correction using COCONUT R package
##  
#########################################################################

## Load data to be co-normalized
load("GSEs_discovery.RData") # follow the data structure mentioned in COCONUT R package

## COCONUT normalization
library(COCONUT)
GSEs.COCONUT <- COCONUT(GSEs=GSEs.discovery,
                        control.0.col="cond",
                        byPlatform=FALSE)
GSEs.COCO.combined <- combineCOCOoutput(GSEs.COCONUT)

