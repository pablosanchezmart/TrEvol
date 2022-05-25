library(testthat)
library(TrEvol)

setwd("tests")

test_check("TrEvol")

### simulate dataset ####

nObs <- 100

tr <- phytools::pbtree(n = nObs)

df <- simulateDataSet(nObs = nObs, phylogeny = tr)

specifications <- defineModelsSpecifications()

initializeTrEvo()

### compute variance partition ####

varianceResults <- computeVariancePartition(traits = c("BM_HC_1", "BM_HC_1") , environmentalVariables = "BM_HC_predictor", dataset = df,
                                            phylogeny = tr, model.specifications = specifications, forceRun = F)


### compute covariance partition ####



### plot results ####
