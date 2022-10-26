# remove(list = ls())
#
# library(testthat)
# library(TrEvol)
#
#
# # test_check("TrEvol")
#
# ### simulate dataset ####
#
# setwd("tests")
#
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
# # specifications <- defineModelsSpecifications(number.iterations = 1000000, burning = 1000, thinning = 10)
# specifications <- defineModelsSpecifications(number.iterations = 100, burning = 10, thinning = 2)
#
# specifications <- defineModelsSpecifications()
#
# initializeTrEvo()
#
# TRAITS <- names(df$data)[!is.na(stringr::str_extract(names(df$data), "trait"))]
# # TRAITS <- TRAITS[1:2]
#
# ENVIRONMENTALVARIABLES <- names(df$data)[!is.na(stringr::str_extract(names(df$data), "env"))]
# ENVIRONMENTALVARIABLES <- ENVIRONMENTALVARIABLES[1]
#
# ### compute variance partition ####
#
# varianceResults_BM_HC_predictor <- computeVariancePartition(traits = TRAITS, environmental.variables = "phylo_G1_envPred", dataset = df$data,
#                                             phylogeny = tr, model.specifications = specifications, force.run = T, save = F)
#
# varianceResults_BM_HC_predictor$varianceResults
#
#
# ### compute covariance partition ####
#
# covarianceResults <- computeCovariancePartition(traits = c(TRAITS, ENVIRONMENTALVARIABLES), dataset = df$data,
#                  phylogeny = tr, model.specifications = specifications, force.run = F, save = F)
#
# covarianceResults_BM_HC_predictor <- computeCovariancePartition(traits = TRAITS, environmental.variables = "phylo_G1_envPred", dataset = df$data,
#                                                                 phylogeny = tr, model.specifications = specifications, force.run = F, save = F)
#
#
# covarianceResults_BM_HC_predictor$covarianceResults
#
# covarianceResults_BM_HC_predictor$covarianceResults$Total_coordination
#
# covarianceResults_BM_HC_predictor$covarianceResults$Pure_coordinated_phylogenetic_conservatism + covarianceResults_BM_HC_predictor$covarianceResults$Coordinated_phylogenetic_niche_conservatism + covarianceResults_BM_HC_predictor$covarianceResults$Pure_environmental_coordination + covarianceResults_BM_HC_predictor$covarianceResults$Residual_coordination
#
# ### compute variance covariance partition ####
#
# varcovResults <- computeVarianceCovariancePartition(traits = TRAITS, dataset = df$data,
#                                                     phylogeny = tr, model.specifications = specifications, force.run = T, save = T, showRelativeResults = T)
#
# varcovResults$covarianceResults
#
# # ### plot results ####
#
# plotData(phylogeny = tr, dataset = df$data, variables = TRAITS)
#
# plotDta <- getPlotNetworkData(covariance.results = covarianceResults_BM_HC_predictor$covarianceResults, variance.results = varianceResults_BM_HC_predictor$varianceResults,
#                               variance.type = "Total_phylogenetic_conservatism", covariance.type = "Total_coordinated_phylogenetic_conservatism", only.significant = T)
# plotDta
#
# ## Two models approach
#
# plotNetwork(covariance.results = covarianceResults_BM_HC_predictor$covarianceResults, variance.results = varianceResults_BM_HC_predictor$varianceResults,
#             variance.type = "Total_phylogenetic_conservatism", covariance.type =  "Total_coordinated_phylogenetic_conservatism",
#             only.significant = T, threshold = 0.1, layout = "circular")
#
# ## One model approach
#
# plotNetwork(covariance.results = varcovResults$covarianceResults, variance.results = varcovResults$varianceResults,
#             variance.type = "Total_phylogenetic_conservatism", covariance.type =  "Total_coordinated_phylogenetic_conservatism",
#             only.significant = T, threshold = 0.1, layout = "circular", displayDegreeAsNodeSize = T)
#
# varcovResults$covarianceResults %>% dplyr::select(Trait_1, Trait_2, Total_coordinated_phylogenetic_conservatism) %>% head()
#
# head(varcovResults$covarianceResults)
#
# ## Compared posterior distributions
#
# plot(varcovResults$individual.models.results$tri_phylo_G1_trait1_phylo_G1_trait2_phylo_G1_envPred$covariancePartitionDistributions$totalCoordinatedPhylogeneticConservatism)
#
# plot(covarianceResults_BM_HC_predictor$individual.models.results$bi_phylo_G1_trait2_phylo_G1_trait1$covariancePartitionDistributions$totalCoordinatedPhylogeneticConservatism)
#
#
# #### IMPUTATION -----------------------------------------------------------------  ####
#
# #### PREDICTIONS --------------------------------------------------------------- ####
#
# # Pablo Sanchez Martinez
#
# dataGroup <- "simulations"
#
# print(paste0("Extracting results for: ", dataGroup))
#
# #### MODELS OUTPUTS ---------------------------------------------------------- ####
#
# for(f in list.files(path = paste0("outputs/outputs_", dataGroup, "/models_outputs/"),  full.names = T)){
#     load(file = f)
# }
#
#
# #### MODELS PERFORMANCE --------------------------------------------------------- ####
#
#
# ### Phylogeny as predictor ####
#
# # tr_pred <- read.tree("data/treevol_phylogeny_genus_lvl_2021.tre")
# df <- df$data
#
# df_pred <- df %>% dplyr::filter(animal %in% tr$tip.label)
# df_pred <- df_pred %>% dplyr::rename(taxon = animal)
#
# tr_pred <- ape::keep.tip(tr, tr$tip.label[tr$tip.label %in% df_pred$taxon])
#
# predictors <- ENVIRONMENTALVARIABLES
#
# for(predictor in predictors){
#     df_pred <- df_pred %>% dplyr::filter(!is.na(predictor))
# }
#
#
# imputationPredictors <- NULL
# propNA <- 0.2
# numberIterations <- 10
# forceRunImputation <- T
#
# TRAITS <- TRAITS[1:2]
#
# # df_pred[1:10, 2] <- NA
# # df_pred[1:10, 3] <- NA
#
# ## Evolutionary order
#
# df_imp <- imputeTraits(imputationVariables = TRAITS[1:2], dataset = df_pred, phylogeny = tr_pred, predictors = ENVIRONMENTALVARIABLES,
#                        orderCriterium = "Total_coordination", prodNAs = propNA,
#                        IterationsNumber = numberIterations, parallelization = T, clustersNumber = 2,
#                        forceRun = forceRunImputation, numberOfImputationRounds = 4)
#
# df_imp$round3$ximp
# head(df_imp$round1$ximp)
#
# df_imp$round3$ximp
#
# errors <- rbind(cbind("type" = "imputation1", df_imp$round1$predictivePerformance[, c("Variable", "NRMSE", "R2")]),
#                 cbind("type" = "imputation2", df_imp$round2$predictivePerformance[, c("Variable", "NRMSE", "R2")]),
#                 cbind("type" = "imputation3", df_imp$round3$predictivePerformance[, c("Variable", "NRMSE", "R2")]),
#                 cbind("type" = "imputation4", df_imp$round4$predictivePerformance[, c("Variable", "NRMSE", "R2")]))
#
# pNRMSE <- ggplot2::ggplot(errors, ggplot2::aes(x = Variable, y = NRMSE, color = type)) +
#   ggplot2::geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
#   ggplot2::geom_boxplot()
# pNRMSE
#
# pR2 <- ggplot2::ggplot(errors, ggplot2::aes(x = Variable, y = R2, color = type)) +
#   ggplot2::geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
#   ggplot2::geom_boxplot()
# pR2
#
#
# ## Random order
#
# df_imp_rand <- imputeTraits(dataset = df_pred, phylogeny = tr_pred, correlationsTraitsResults = NULL,
#                        varianceResults = NULL,
#                        orderCriterium = NULL,
#                        imputationVariables = variablesToImpute, predictors = ENVIRONMENTALVARIABLES, numberOfPhyloCoordinates = 5, prodNAs = propNA,
#                        IterationsNumber = numberIterations, clustersNumber = 2,
#                        forceRun = forceRunImputation)
#
# df_imp_rand$predictivePerformance2
#
# errors_rand <- rbind(cbind("type" = "imputation1", df_imp_rand$predictivePerformance_all_iterations[, c("Variable", "NRMSE", "R2")]),
#                 cbind("type" = "imputation2", df_imp_rand$predictivePerformance_all_iterations2[, c("Variable", "NRMSE", "R2")]))
#
# pNRMSE_rand <- ggplot2::ggplot(errors_rand, ggplot2::aes(x = Variable, y = NRMSE, color = type)) +
#   ggplot2::geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
#   ggplot2::geom_boxplot()
# pNRMSE_rand
#
# pR2_rand <- ggplot2::ggplot(errors_rand, ggplot2::aes(x = Variable, y = R2, color = type)) +
#   ggplot2::geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
#   ggplot2::geom_boxplot()
# pR2_rand
#
#
# df_imp$predictivePerformance2$R2 - df_imp_rand$predictivePerformance$R2
#
# df_imp$predictivePerformance$R2 - df_imp_rand$predictivePerformance$R2
#
# #
# # # Delete files
#
# unlink("outputs", recursive = T)
# unlink("results", recursive = T)
