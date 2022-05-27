# library(testthat)
# library(TrEvol)

setwd("tests")

# test_check("TrEvol")

### simulate dataset ####

# devtools::load_all()
#
# nObs <- 100
#
# tr <- phytools::pbtree(n = nObs)
#
# df <- simulateDataSet(nObs = nObs, phylogeny = tr, vcvMatrix = matrix(c(1, 0.9, 0.8, 0, 0.1, 0.2,
#                                                                         0.9, 1, 0.8, 0, 0.1, 0.2,
#                                                                         0.8, 0.8, 1, 0, 0.1, 0.2,
#                                                                         0, 0, 0, 1, -0.9, -0.8,
#                                                                         0.1, 0.1, 0.1, -0.9, 1, 0.8,
#                                                                         0.2, 0.2, 0.2, -0.8, 0.8, 1), ncol = 6))
#
# specifications <- defineModelsSpecifications(n_itterations = 100000, burning = 100, thinning = 5)
#
# initializeTrEvo()
#
# TRAITS <- c("BM_HC_1", "BM_HC_2", "BM_LC_1", "BM_LC_2", "nonBM_HC_1", "nonBM_HC_2", "nonBM_LC_1", "nonBM_LC_2")
#
# ENVIRONMENTALVARIABLES <- c("BM_HC_predictor", "BM_LC_predictor", "nonBM_HC_predictor", "nonBM_LC_predictor")
#
# ### compute variance partition ####
#
# varianceResults_BM_HC_predictor <- computeVariancePartition(traits = TRAITS, environmentalVariables = "BM_HC_predictor", dataset = df,
#                                             phylogeny = tr, model.specifications = specifications, forceRun = T)
#
# varianceResults_nonBM_HC_predictor <- computeVariancePartition(traits = TRAITS , environmentalVariables = "nonBM_HC_predictor", dataset = df,
#                                             phylogeny = tr, model.specifications = specifications, forceRun = T)
#
# varianceResults_BM_HC_predictor$varianceResults
# varianceResults_nonBM_HC_predictor$varianceResults
#
# varianceResults$varianceResults$Pure_phylogenetic_conservatism + varianceResults$varianceResults$Phylogenetic_niche_conservatism + varianceResults$varianceResults$Pure_environmental + varianceResults$varianceResults$Residual
#
# ### compute covariance partition ####
#
# covarianceResults_BM_HC_predictor <- computeCovariancePartition(traits = TRAITS, environmentalVariables = "BM_HC_predictor", dataset = df,
#                                                                 phylogeny = tr, model.specifications = specifications, forceRun = F)
#
# covarianceResults_nonBM_HC_predictor <- computeCovariancePartition(traits = TRAITS, environmentalVariables = "nonBM_HC_predictor", dataset = df,
#                                                                 phylogeny = tr, model.specifications = specifications, forceRun = F)
#
#
# covarianceResults_BM_HC_predictor$covarianceResults
#
# covarianceResults_BM_HC_predictor$covarianceResults$Total_coordination
#
# covarianceResults_BM_HC_predictor$covarianceResults$Pure_coordinated_phylogenetic_conservatism + covarianceResults_BM_HC_predictor$covarianceResults$Coordinated_niche_phylogenetic_conservatism + covarianceResults_BM_HC_predictor$covarianceResults$Pure_environmental_coordination + covarianceResults_BM_HC_predictor$covarianceResults$Residual_coordination
#
# ### plot results ####
#

# unlink("outputs", recursive = T)
# unlink("results", recursive = T)
