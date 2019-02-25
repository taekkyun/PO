## Load data to be co-normalized
load("GSEs_discovery.RData")

## COCONUT normalization
library(COCONUT)
GSEs.COCONUT <- COCONUT(GSEs=GSEs.discovery,
                        control.0.col="cond",
                        byPlatform=FALSE)
GSEs.COCO.combined <- combineCOCOoutput(GSEs.COCONUT)

