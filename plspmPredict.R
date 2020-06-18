##########################################################################################
# plspmPredict.R
# Description: This function predicts PLS-PM latent and measurement variables from
#              a 'plspm' object ('plspm' R-package)
# This is based on the publication:
#   Shmueli, G., Ray, S., Estrada, J. M., & Chatla, S. (n.d.). The Elephant in the Room:
#   Evaluating the Predictive Performance of Partial Least Squares (PLS) Path Models (2015).
#   SSRN Electronic Journal SSRN Journal.
#
#  The scriptis based on the following code:
#   https://github.com/ISS-Analytics/pls-predict/blob/master/lib/PLSpredict.R
#   in order to work directly with an plspm object and to allow for raster data predicitons
#
# data input:
# -----------
#   - pls: An plspm object from the plspm package
#   - dat: dataframe or Raster Stack with the model predictors
#
# @author: Javier Lopatin
#
##########################################################################################

plspmPredict <- function(pls, dat)
{
  if (class(pls) != "plspm")
    stop("\n'plspmPredict()' requires a 'plspm' object")
  if (all(class(dat) != "data.frame" && class(dat) != "RasterStack" && class(dat) != "RasterBrick"))
    stop("\n'plspmPredict()' requires a 'data.frame', 'RasterStack', or 'RasterBrick' object")
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
  path_coef   <- pls$path_coefs # path coefficients

  # Extract and Normalize the measurements for the model
  normDataTrain <- scale(pls$data[, mmVariables], TRUE, TRUE)
  # Extract Mean and Standard Deviation of measurements for future prediction
  meanData <- attr(normDataTrain, "scaled:center")
  sdData   <- attr(normDataTrain, "scaled:scale")

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
    exVar = strsplit(as.character(pls$effects[,1][i]), "[->]")[[1]][1]
    enVar = strsplit(as.character(pls$effects[,1][i]), "[->]")[[1]][3]
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

  # =======================================================
  # predict PLS-PM values
  # =======================================================

  # Extract Measurements needed for Predictions
  if (method == 'dat'){
    normData <- dat[, pMeasurements]
    # Normalize data
    for (i in pMeasurements){
      normData[,i] <- (normData[,i] - meanData[i])/sdData[i]
    }
    # Convert dataset to matrix
    normData <- data.matrix(normData)

    # Add empty columns to normData for the estimated measurements
    for (i in 1:length(eMeasurements)){
      normData = cbind(normData, seq(0,0,length.out = nrow(normData)))
      colnames(normData)[length(colnames(normData))] = eMeasurements[i]
    }

    # Estimate Factor Scores from Outer Path
    fscores <- normData %*% outer_weights
    # Estimate Factor Scores from Inner Path and complete Matrix
    fscores <- fscores + fscores %*% t(path_coef)

    # Predict Measurements with loadings
    predictedMeasurements <- fscores %*% t(outer_loadings)

    # Denormalize data
    for (i in mmVariables){
      predictedMeasurements[,i]<-(predictedMeasurements[,i] * sdData[i]) + meanData[i]
    }

    # Calculating the measurement residuals and accuracies if validation data is provided
    if (!is.na(sum( match( eMeasurements, colnames(dat) ) ))){ # if validation data is presented for the endogenous variables

      # measurement variables data
      mmData = dat[, eMeasurements]

      # get residuals
      mmResiduals <- dat[,resMeasurements] - predictedMeasurements[,resMeasurements]

      # get accuracies of measurement predictions
      fit_measurements <- matrix(ncol=length(resMeasurements), nrow = 4, byrow = T)
      colnames(fit_measurements) <- resMeasurements
      rownames(fit_measurements) <- c('r_square','RMSE','nRMSE','bias')
      for (i in 1:length(resMeasurements)){
        fit_measurements[1,i] <- cor(dat[,resMeasurements[i]], predictedMeasurements[,resMeasurements[i]], method="pearson")
        fit_measurements[2,i] <- sqrt(mean((dat[,resMeasurements[i]] - predictedMeasurements[,resMeasurements[i]])^2))
        fit_measurements[3,i] <- (fit_measurements[2,i]/(max(dat[,resMeasurements[i]]) - min(dat[,resMeasurements[i]])))*100
        lm = lm(predictedMeasurements[,resMeasurements[i]] ~ dat[,resMeasurements[i]]-1)
        fit_measurements[4,i] <- 1-coef(lm)
      }
    } else {
      mmData = dat[,pMeasurements]
      residuals = NA
      fit_measurements = NA
    }

    # Prepare return Object
    predictResults <- list(mmData = mmData,
                           mmPredicted = predictedMeasurements[,resMeasurements],
                           mmResiduals = mmResiduals,
                           mmfit = fit_measurements,
                           Scores = fscores)
  }

  if (method == 'rst'){
    names(dat) = pMeasurements
    normData <- (dat[[1]] - meanData[1])/sdData[1]
    for (i in 2:length(pMeasurements)){
      #normData[[i]] <- (dat[[i]] - meanData[i])/sdData[i]
      normData <- addLayer(normData, (dat[[i]] - meanData[i])/sdData[i])
    }

    # Add empty columns to normData for the estimated measurements
    for (i in 1:length(eMeasurements)){
      normData <- addLayer(normData, normData[[1]])
      names(normData)[length(names(normData))] <- eMeasurements[i]
      normData[[length(names(normData))]][normData[[length(names(normData))]]] <- 0
    }

    # Estimate Factor Scores from Outer Path
    fscores <- raster(normData)
    for (i in 1:length(ltVariables)){
      fscores <- addLayer(fscores, calc(normData[[which(outer_weights[,i] == 1)]], fun=prod))
    }
    # Estimate Factor Scores from Inner Path and complete Matrix
    for (i in 1:length(ltVariables)){
      idx = which(t(path_coef)[,i] != 0)
      n = length(idx)
      if (n == 0){
        fscores[[i]] <- fscores[[i]]
        names(fscores[[i]]) = ltVariables[i]
      } else if (n == 1){
        fscores[[i]] <- fscores[[i]] + (fscores[[idx]] * t(path_coef)[,i][idx])
        names(fscores[[i]]) = ltVariables[i]
      } else {
        r = raster(fscores)
        for (j in 1:n){ r = addLayer(r, (fscores[[idx[j]]] * t(path_coef)[,i][idx][j])) }
        fscores[[i]] <- fscores[[i]] + calc(r, fun=sum)
        names(fscores[[i]]) = ltVariables[i]
      }
    }

    # Predict Measurements with loadings
    predictedMeasurements <- raster(fscores[[1]])

    if (length(eMeasurements) == 1){
      idx = which(mmVariables == eMeasurements) # measurement
      idx2 = which(ltVariables == enVariables) # LV
      predictedMeasurements <- addLayer(predictedMeasurements, fscores[[idx2]] * t(outer_loadings)[idx2, idx])
      names(predictedMeasurements) = eMeasurements[i]
    } else{
      idx2 = which(ltVariables == enVariables) # LV
      for (i in 1:length(eMeasurements)){ # measurements
        idx = which(mmVariables == eMeasurements[i]) # measurement
        r = raster(fscores[[1]])
        for (k in 1:length(idx)){ r = addLayer(r, (fscores[[idx2]] * t(outer_loadings)[idx2, idx])) }
        predictedMeasurements <- addLayer(predictedMeasurements, calc(r, fun=sum))
        names(predictedMeasurements)[i] = eMeasurements[i]
     }
    }

    # Denormalize data
    for (i in 1:length(eMeasurements)){
      idx = which(names(sdData) == eMeasurements[i])
      predictedMeasurements[[i]] <- (predictedMeasurements[[i]] * sdData[idx]) + meanData[idx]
    }

    # Prepare return Object
    predictResults <- list(mmPredicted = predictedMeasurements,
                           Scores = fscores)

  }

  # =======================================================
  # function output
  # =======================================================

  class(predictResults) <- "plspmPredict"
  return(predictResults)

}
