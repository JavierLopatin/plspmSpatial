##########################################################################################
# groupsPredict.R
# Description: This function predicts PLS-PM latent and measurement variables from
#              a 'plspm.groups' object ('plspm' R-package)
# This is based on the publication:
#   Shmueli, G., Ray, S., Estrada, J. M., & Chatla, S. (n.d.). The Elephant in the Room:
#   Evaluating the Predictive Performance of Partial Least Squares (PLS) Path Models (2015).
#   SSRN Electronic Journal SSRN Journal.
#
#
# data input:
# -----------
#   - pls: An plspm object from the plspm package
#   - pls.groups: An plspm.groups object from the plspm package
#   - train.groups: A vector with the classes of the dataset used to train plspm
#   - dat: dataframe or Raster Stack with the model predictors
#
# @author: Javier Lopatin
#
##########################################################################################


plspm.groupsPredict <- function(pls, pls.groups, train.groups, dat){
  if (class(pls) != "plspm")
    stop("\n'plspm.groupsPredict()' requires a 'plspm' object")
  if (all(class(dat) != "data.frame" && class(dat) != "RasterStack" && class(dat) != "RasterBrick"))
    stop("\n'plspm.groupsPredict()' requires a 'data.frame', 'RasterStack', or 'RasterBrick' object")
  # Determine method
  if (class (dat) ==  "data.frame") method <- 'dat'
  if (class (dat) == 'RasterBrick') method <- 'rst'
  if (class (dat) == 'RasterStack') method <- 'rst'


  # =======================================================
  # inputs setting
  # =======================================================

  # from plspm object
  ltVariables <- pls$model$gen$lvs_names # latent variables
  mmVariables <- pls$model$gen$mvs_names # measurement variables
  group1_obs <- which(train.groups == unique(train.groups)[1])
  group2_obs <- which(train.groups == unique(train.groups)[2])

  # Extract and Normalize the measurements for the model
  a <- scale(pls$data[, mmVariables], TRUE, TRUE)
  b <- scale(pls$data[group1_obs, mmVariables], TRUE, TRUE)
  c <- scale(pls$data[group2_obs, mmVariables], TRUE, TRUE)
  # Extract Mean and Standard Deviation of measurements for future prediction
  meanData_global <- attr(a, "scaled:center")
  sdData_global   <- attr(a, "scaled:scale")
  meanData_group1 <- attr(b, "scaled:center")
  sdData_group1   <- attr(b, "scaled:scale")
  meanData_group2 <- attr(c, "scaled:center")
  sdData_group2   <- attr(c, "scaled:scale")

  # =======================================================
  # prepare data
  # =======================================================

  # get relationship matrix
  mmMatrix = matrix(nrow = length(mmVariables),
                    ncol = 2, byrow =TRUE,
                    dimnames = list(1:length(mmVariables), c("latent","measurement")))
  mmMatrix[,'latent'] = as.character(pls$outer_model[,2])
  mmMatrix[,'measurement'] = as.character(pls$outer_model[,1])

  # Create a matrix of outer_weights
  outer_weights <- matrix(data=0,
                          nrow=length(mmVariables),
                          ncol=length(ltVariables),
                          dimnames = list(mmVariables,ltVariables))

  #Initialize outer_weights matrix with value 1 for each relationship in the measurement model
  for (i in 1:length(ltVariables))  {
    outer_weights [mmMatrix[mmMatrix[,"latent"]==ltVariables[i],
                            "measurement"],
                   ltVariables[i]] = 1
  }

  # get relationship matrix
  smMatrix = matrix(nrow = length(pls$effects[,1]),
                    ncol = 2, byrow =TRUE,
                    dimnames = list(1:length(pls$effects[,1]), c("source","target")))

  for (i in 1:length(pls$effects[,1])){
    exVar = strsplit(as.character(pls$effects[,1][i]), "[-->]")[[1]][1]
    enVar = strsplit(as.character(pls$effects[,1][i]), "[-->]")[[1]][3]
    smMatrix[i,1] <- gsub(" ", "", exVar, fixed = TRUE)
    smMatrix[i,2] <- gsub(" ", "", enVar, fixed = TRUE)
  }

  # Create a matrix of outer_loadins
  outer_loadings <- matrix(data=0,
                           nrow=length(mmVariables),
                           ncol=length(ltVariables),
                           dimnames = list(mmVariables,ltVariables))

  for (i in 1:length(ltVariables))  {
    mesVar = mmMatrix[mmMatrix[,"latent"]==ltVariables[i], "measurement"]
    idx = which(pls$outer_model$name %in% mesVar)
    outer_loadings[mesVar , ltVariables[i]] = pls$outer_model$loading[idx]
  }

  # Identify Exogenous and Endogenous Variables
  exVariables <- unique(smMatrix[,1])
  pMeasurements <- NULL
  for (i in 1:length(exVariables)){
    pMeasurements <- c(pMeasurements,mmMatrix[mmMatrix[,"latent"]==exVariables[i],"measurement"])
  }

  enVariables <- unique(smMatrix[,2])
  resMeasurements <- NULL
  for (i in 1:length(enVariables)){
    resMeasurements <- c(resMeasurements, mmMatrix[mmMatrix[, "latent"] == enVariables[i],"measurement"])
  }

  enVariables <- setdiff(enVariables,exVariables)
  eMeasurements <- NULL
  for (i in 1:length(enVariables)){
    eMeasurements <- c(eMeasurements, mmMatrix[mmMatrix[, "latent"] == enVariables[i],"measurement"])
  }

  # Get path coefficients from pls.groups$test
  path_coef_global = matrix(0, nrow=length(ltVariables),ncol=length(ltVariables),
                            dimnames = list(1:length(ltVariables), ltVariables))
  rownames(path_coef_global) <- ltVariables

  path_coef_group1 = matrix(0, nrow=length(ltVariables),ncol=length(ltVariables),
                            dimnames = list(1:length(ltVariables), ltVariables))
  rownames(path_coef_group1) <- ltVariables

  path_coef_group2 = matrix(0, nrow=length(ltVariables),ncol=length(ltVariables),
                            dimnames = list(1:length(ltVariables), ltVariables))
  rownames(path_coef_group2) <- ltVariables

  # global
  for (i in 1:length(ltVariables)){ # iter through the source variables
    for (j in 1:length(ltVariables)){ # iter through the target variables
      exVar = ltVariables[i]
      enVar = ltVariables[j]
      #while (a != exVar && b != enVar){}
      for (k in 1:nrow(pls.groups$test)){ # iter through the alternative of table effects
        a = strsplit(as.character(rownames(pls.groups$test)[k]), "[-->]")[[1]][1]
        b = strsplit(as.character(rownames(pls.groups$test)[k]), "[-->]")[[1]][3]
        if (a == exVar && b == enVar){
          path_coef_global[j,i] <- pls.groups$test$global[k]
        }
      }
    }
  }

  # group1
  for (i in 1:length(ltVariables)){ # iter through the source variables
    for (j in 1:length(ltVariables)){ # iter through the target variables
      exVar = ltVariables[i]
      enVar = ltVariables[j]
      #while (a != exVar && b != enVar){}
      for (k in 1:nrow(pls.groups$test)){ # iter through the alternative of table effects
        a = strsplit(as.character(rownames(pls.groups$test)[k]), "[-->]")[[1]][1]
        b = strsplit(as.character(rownames(pls.groups$test)[k]), "[-->]")[[1]][3]
        if (a == exVar && b == enVar){
          path_coef_group1[j,i] <- pls.groups$test[k,2]
        }
      }
    }
  }

  # group2
  for (i in 1:length(ltVariables)){ # iter through the source variables
    for (j in 1:length(ltVariables)){ # iter through the target variables
      exVar = ltVariables[i]
      enVar = ltVariables[j]
      #while (a != exVar && b != enVar){}
      for (k in 1:nrow(pls.groups$test)){ # iter through the alternative of table effects
        a = strsplit(as.character(rownames(pls.groups$test)[k]), "[-->]")[[1]][1]
        b = strsplit(as.character(rownames(pls.groups$test)[k]), "[-->]")[[1]][3]
        if (a == exVar && b == enVar){
          path_coef_group2[j,i] <- pls.groups$test[k,3]
        }
      }
    }
  }

  # =======================================================
  # predict PLS-PM values
  # =======================================================


  # Extract Measurements needed for Predictions
  if (method == 'dat'){ normData <- dat[, pMeasurements] }

  # Normalize data
  for (i in pMeasurements){
    normData[,i] <-(normData[,i] - meanData[i])/sdData[i]
  }

  # Convert dataset to matrix
  normData <- data.matrix(normData)

  # Add empty columns to normData for the estimated measurements
  for (i in 1:length(eMeasurements)){
    normData = cbind(normData, seq(0,0,length.out =nrow(normData)))
    colnames(normData)[length(colnames(normData))]=eMeasurements[i]
  }

  # Estimate Factor Scores from Outer Path
  fscores <- normData %*% outer_weights

  # Estimate Factor Scores from Inner Path and complete Matrix
  fscores_global <- fscores + fscores %*% t(path_coef_global)
  fscores_group1 <- fscores + fscores %*% t(path_coef_group1)
  fscores_group2 <- fscores + fscores %*% t(path_coef_group2)

  # make a list with fscores
  Scores <- list(globa=fscores_global,
                    group1=fscores_group1,
                    group2=fscores_group2)

  # Predict Measurements with loadings
  predic_global <- fscores_global %*% t(outer_loadings)
  predic_group1 <- fscores_group1 %*% t(outer_loadings)
  predic_group2 <- fscores_group2 %*% t(outer_loadings)

  # Denormalize data
  for (i in mmVariables){
    predic_global[,i] <- (predic_global[,i] * sdData_global[i]) + meanData_global[i]
    predic_group1[,i] <- (predic_group1[,i] * sdData_group1[i]) + meanData_group1[i]
    predic_group2[,i] <- (predic_group2[,i] * sdData_group2[i]) + meanData_group2[i]
  }

  # make a list with predictions
  predictions <- list(globa=predic_global,
                      group1=predic_group1,
                      group2=predic_group2)

  # Calculating the measurement residuals and accuracies if validation data is provided
  if (!is.na(sum( match( eMeasurements, colnames(dat) ) ))){

      # measurement variables data
      mmData = dat[, eMeasurements]

      # get residuals
      mmResiduals_global <- dat[,resMeasurements] - predic_global[,resMeasurements]
      mmResiduals_group1 <- dat[, resMeasurements] - predic_group1[,resMeasurements]
      mmResiduals_group2 <- dat[, resMeasurements] - predic_group2[,resMeasurements]

      # make a list with residuals
      residuals <- list(globa=mmResiduals_global,
                        group1=mmResiduals_group1,
                        group2=mmResiduals_group2)

      # get r-sqaure values
      r_square <- matrix(ncol = length(resMeasurements), nrow = 3, byrow = T)
      colnames(r_square) <- resMeasurements
      rownames(r_square) <- c('global','group1','group2')
      for (i in 1:length(resMeasurements)){
        r_square[1,i] <- cor(dat[,resMeasurements[i]], predic_global[,resMeasurements[i]], method="pearson")
        r_square[2,i] <- cor(dat[,resMeasurements[i]], predic_group1[,resMeasurements[i]], method="pearson")
        r_square[3,i] <- cor(dat[,resMeasurements[i]], predic_group2[,resMeasurements[i]], method="pearson")
      }

      # get RMSE values
      RMSE <- matrix(ncol = length(resMeasurements), nrow = 3, byrow = T)
      colnames(RMSE) <- resMeasurements
      rownames(RMSE) <- c('global','group1','group2')
      for (i in 1:length(resMeasurements)){
        RMSE[1,i] <- sqrt(mean((dat[,resMeasurements[i]] - predic_global[,resMeasurements[i]])^2))
        RMSE[2,i] <- sqrt(mean((dat[,resMeasurements[i]] - predic_group1[,resMeasurements[i]])^2))
        RMSE[3,i] <- sqrt(mean((dat[,resMeasurements[i]] - predic_group2[,resMeasurements[i]])^2))
      }

      # get normalize-RMSE values
      nRMSE <- matrix(ncol = length(resMeasurements), nrow = 3, byrow = T)
      colnames(nRMSE) <- resMeasurements
      rownames(nRMSE) <- c('global','group1','group2')
      for (i in 1:length(resMeasurements)){
        nRMSE[1,i] <- (RMSE[1,i]/( max(dat[,resMeasurements[i]]) - min(dat[,resMeasurements[i]]) ))*100
        nRMSE[2,i] <- (RMSE[2,i]/( max(dat[,resMeasurements[i]]) - min(dat[,resMeasurements[i]]) ))*100
        nRMSE[3,i] <- (RMSE[3,i]/( max(dat[,resMeasurements[i]]) - min(dat[,resMeasurements[i]]) ))*100
      }

      # get bias values
      bias <- matrix(ncol = length(resMeasurements), nrow = 3, byrow = T)
      colnames(bias) <- resMeasurements
      rownames(bias) <- c('global','group1','group2')
      for (i in 1:length(resMeasurements)){
        lm1 = lm(predic_global[,resMeasurements[i]] ~ dat[,resMeasurements[i]]-1)
        lm2 = lm(predic_group1[,resMeasurements[i]] ~ dat[,resMeasurements[i]]-1)
        lm3 = lm(predic_group2[,resMeasurements[i]] ~ dat[,resMeasurements[i]]-1)
        bias[1,i] <- 1-coef(lm1)
        bias[2,i] <- 1-coef(lm2)
        bias[3,i] <- 1-coef(lm3)
      }
  } else {
    mmData = dat[,pMeasurements]
    residuals = NA
    r_square = NA
    RMSE = NA
    nRMSE = NA
    bias = NA
  }

  # Prepare return Object
  predictResults <- list(mmData = dat[,resMeasurements],
                         mmPredicted = predictions,
                         mmResiduals = residuals,
                         Scores = Scores,
                         r_square = r_square,
                         RMSE = RMSE,
                         nRMSE = nRMSE,
                         bias = bias)

  class(predictResults) <- "plspmPredict"
  return(predictResults)

}
