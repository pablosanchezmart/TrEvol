randomForestImpute <- function(xmis, maxiter = 10, ntree = 100,
                       mtry = floor(sqrt(ncol(xmis))), replace = TRUE,
                       classwt = NULL, cutoff = NULL, strata = NULL,
                       sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                       xtrue = NA, parallelize = FALSE)
{
  require(missForest)
  require(doRNG)
  require(randomForest)
  combine <- randomForest::combine
  ## stop in case of wrong inputs passed to randomForest
  n <- nrow(xmis)
  # p <- ncol(xmis)

  ## extract missingness pattern

  NAloc <- is.na(xmis)        # where are missings
  noNAvar <- apply(NAloc, 2, sum) # how many are missing in the vars
  varsToImpute <- names(which(noNAvar != 0))

  p <- length(varsToImpute)

  if (!is.null(classwt)){
    stopifnot(length(classwt) == p, typeof(classwt) == 'list')
  }
  if (!is.null(cutoff)){
    stopifnot(length(cutoff) == p, typeof(cutoff) == 'list')
  }
  if (!is.null(strata)){
    stopifnot(length(strata) == p, typeof(strata) == 'list')
  }
  if (!is.null(nodesize)){
    stopifnot(length(nodesize) == 2)
  }

  ## remove completely missing variables
  if (any(apply(is.na(xmis), 2, sum) == n)){
    indCmis <- which(apply(is.na(xmis), 2, sum) == n)
    xmis <- xmis[,-indCmis]
    # p <- ncol(xmis)
    cat('  removed variable(s)', indCmis,
        'due to the missingness of all entries\n')
  }

  ## return feedback on parallelization setup
  if (parallelize) {
    if (foreach::getDoParWorkers() == 1) {
      stop("You must register a 'foreach' parallel backend to run 'missForest' in parallel. Set 'parallelize' to 'no' to compute serially.")
    }
    # if (foreach::getDoParWorkers() > p){
    #   stop("The number of parallel cores should not exceed the number of variables (p=", p, ")")
    # }
  }

  ## perform initial S.W.A.G. on xmis (mean imputation)
  ximp <- xmis
  varType <- character(p)
  names(varType) <- varsToImpute
  for (t.co in 1:p) {
    if (is.numeric(xmis[[t.co]])) {
      varType[t.co] <- 'numeric'
      next()
    }
    if (is.factor(xmis[[t.co]])) {
      varType[t.co] <- 'factor'
      next()
    }
    stop(sprintf('column %s must be factor or numeric, is %s', names(xmis)[t.co], class(xmis[[t.co]])))
  }

  if(length(varsToImpute) < 1){
    stop("No missing values in dataset")
  }
  sort.j <- c(1:length(noNAvar))

  sort.noNAvar <- noNAvar[sort.j]

  ## compute a list of column indices for variable parallelization
  nzsort.j <- sort.j[sort.noNAvar > 0]

  ## force column loop to be sequential
  '%cols%' <- get('%do%')
  idxList <- nzsort.j

  ## output
  Ximp <- vector('list', maxiter)

  ## initialize parameters of interest
  iter <- 1
  k <- length(unique(varType))
  convNew <- rep(0, k)
  convOld <- rep(Inf, k)

  OOBerror <- rep(NA, length(varsToImpute))
  names(OOBerror) <- varsToImpute

  ## setup convergence variables w.r.t. variable types
  if (k == 1){
    if (unique(varType) == 'numeric'){
      names(convNew) <- c('numeric')
    } else {
      names(convNew) <- c('factor')
    }
    convergence <- c()
    OOBerr <- numeric(1)
  } else {
    names(convNew) <- c('numeric', 'factor')
    convergence <- matrix(NA, ncol = 2)
    OOBerr <- numeric(2)
  }

  ## function to yield the stopping criterion in the following 'while' loop
  stopCriterion <- function(varType, convNew, convOld, iter, maxiter){
    k <- length(unique(varType))
    if (k == 1){
      (convNew < convOld) & (iter < maxiter)
    } else {
      ((convNew[1] < convOld[1]) | (convNew[2] < convOld[2])) & (iter < maxiter)
    }
  }

  ximp.old <- ximp

  ## iterate missForest
  while (stopCriterion(varType, convNew, convOld, iter, maxiter)){
    # print(iter)
    if (iter > 2){
      convOld <- convNew
      OOBerrOld <- OOBerr
    }

    t.start <- proc.time()

    ximp <- ximp[, -which(colnames(ximp) %in% varsToImpute)]

    for (varInd in varsToImpute) {

      ## Select only variables without missing values

      ximp <- cbind(ximp.old[, varInd], ximp)
      colnames(ximp)[1] <- varInd


      if (noNAvar[[varInd]] != 0) {
          obsi <- !NAloc[, varInd]
          misi <- NAloc[, varInd]
          obsY <- as.matrix(ximp[obsi, varInd])
          obsX <- as.matrix(ximp[obsi, -which(colnames(ximp) == varInd)])
          misX <- as.matrix(ximp[misi, -which(colnames(ximp) == varInd)])
          typeY <- varType[varInd]
          if (typeY == "numeric") {
            if (parallelize) {
              xntree <- NULL
              RF <- foreach::foreach(xntree = iterators:::idiv(ntree, chunks = foreach::getDoParWorkers()),
                            .combine = "combine", .multicombine = TRUE,
                            .packages = 'randomForest') %dorng% {
                              randomForest( x = obsX,
                                            y = obsY,
                                            ntree = xntree,
                                            mtry = mtry,
                                            replace = replace,
                                            sampsize = if (!is.null(sampsize)) sampsize[[varInd]] else
                                              if (replace) nrow(obsX) else ceiling(0.632 * nrow(obsX)),
                                            nodesize = if (!is.null(nodesize)) nodesize[1] else 1,
                                            maxnodes = if (!is.null(maxnodes)) maxnodes else NULL)
                            }
              ## record out-of-bag error
              OOBerror[varInd] <- mean((predict(RF) - RF$y) ^ 2, na.rm = TRUE)
              # OOBerror[varInd] <- RF$mse[ntree]
            } else {
              RF <- randomForest( x = obsX,
                                  y = obsY,
                                  ntree = ntree,
                                  mtry = mtry,
                                  replace = replace,
                                  sampsize = if (!is.null(sampsize)) sampsize[[varInd]] else
                                    if (replace) nrow(obsX) else ceiling(0.632 * nrow(obsX)),
                                  nodesize = if (!is.null(nodesize)) nodesize[1] else 1,
                                  maxnodes = if (!is.null(maxnodes)) maxnodes else NULL)
              ## record out-of-bag error
              OOBerror[varInd] <- RF$mse[ntree]
            }
            misY <- predict(RF, misX)
          } else {
            obsY <- factor(obsY)
            summarY <- summary(obsY)
            if (length(summarY) == 1) {
              misY <- factor(rep(names(summarY), sum(misi)))
            } else {
              if (parallelize) {
                RF <- foreach::foreach(xntree = iterators:::idiv(ntree, chunks = foreach::getDoParWorkers()),
                              .combine = "combine", .multicombine = TRUE,
                              .packages = 'randomForest') %dorng% {
                                randomForest(
                                  x = obsX,
                                  y = obsY,
                                  ntree = xntree,
                                  mtry = mtry,
                                  replace = replace,
                                  classwt = if (!is.null(classwt)) classwt[[varInd]] else
                                    rep(1, nlevels(obsY)),
                                  cutoff = if (!is.null(cutoff)) cutoff[[varInd]] else
                                    rep(1/nlevels(obsY), nlevels(obsY)),
                                  strata = if (!is.null(strata)) strata[[varInd]] else obsY,
                                  sampsize = if (!is.null(sampsize)) sampsize[[varInd]] else
                                    if (replace) nrow(obsX) else ceiling(0.632 * nrow(obsX)),
                                  nodesize = if (!is.null(nodesize)) nodesize[2] else 5,
                                  maxnodes = if (!is.null(maxnodes)) maxnodes else NULL)
                              }
                ## record out-of-bag error
                ne <- as.integer(predict(RF)) != as.integer(RF$y)
                ne <- ne[! is.na(ne)]
                OOBerror[varInd] <- sum(ne) / length(ne)
              } else {
                RF <- randomForest::randomForest(x = obsX,
                                   y = obsY,
                                   ntree = ntree,
                                   mtry = mtry,
                                   replace = replace,
                                   classwt = if (!is.null(classwt)) classwt[[varInd]] else
                                     rep(1, nlevels(obsY)),
                                   cutoff = if (!is.null(cutoff)) cutoff[[varInd]] else
                                     rep(1 / nlevels(obsY), nlevels(obsY)),
                                   strata = if (!is.null(strata)) strata[[varInd]] else obsY,
                                   sampsize = if (!is.null(sampsize)) sampsize[[varInd]] else
                                     if (replace) nrow(obsX) else ceiling(0.632 * nrow(obsX)),
                                   nodesize = if (!is.null(nodesize)) nodesize[2] else 5,
                                   maxnodes = if (!is.null(maxnodes)) maxnodes else NULL)
                ## record out-of-bag error
                OOBerror[varInd] <- RF$err.rate[[ntree, 1]]
              }
              ## predict missing parts of Y
              misY <- predict(RF, misX)
            }
          }
          ximp[misi, varInd] <- misY
        }
      }

    Ximp[[iter]] <- ximp

    t.co2 <- 1
    ## check the difference between iteration steps
    if(iter > 1){

      for (t.type in names(convNew)){
        t.ind <- which(varType == t.type)
        if (t.type == 'numeric'){
          convNew[t.co2] <- sum((ximp[, t.ind] -  Ximp[[iter-1]][, t.ind] )^2) / sum(ximp[, t.ind]^2)
        } else {
          dist <- sum(as.character(as.matrix(ximp[, t.ind])) != as.character(as.matrix( Ximp[[iter-1]][, t.ind])))
          convNew[t.co2] <- dist / (n * sum(varType == 'factor'))
        }
        t.co2 <- t.co2 + 1
      }

    }


    ## compute estimated imputation error
      OOBerr <- OOBerror
      varType <- varType[varsToImpute]
      names(OOBerr)[varType == 'numeric'] <- 'MSE'
      names(OOBerr)[varType == 'factor'] <- 'PFC'

    if (any(!is.na(xtrue))){
      err <- suppressWarnings(mixError(ximp, xmis, xtrue))
    }
      iter <- iter + 1

  }#end while((convNew<convOld)&(iter<maxiter)){

  ## produce output w.r.t. stopping rule
  if (iter == maxiter){
    if (any(is.na(xtrue))){
      out <- list(ximp = Ximp[[iter]], OOBerror = OOBerrOld)
    } else {
      out <- list(ximp = Ximp[[iter]], OOBerror = OOBerr, error = err)
    }
  } else {
    if (any(is.na(xtrue))){
      out <- list(ximp = Ximp[[iter - 1]], OOBerror = OOBerrOld)
    } else {
      out <- list(ximp = Ximp[[iter - 1]], OOBerror = OOBerrOld,
                  error = suppressWarnings(mixError(Ximp[[iter - 1]], xmis, xtrue)))
    }
  }
  class(out) <- 'missForest'
  return(out)
}
