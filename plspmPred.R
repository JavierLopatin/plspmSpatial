##########################################################################################
# plspmPred.R
# Description: This function predicts PLS-PM latent and measurement variables from
#              a 'plspm' object ('plspm' R-package)
# This is based on the publication:
#   Shmueli, G., Ray, S., Estrada, J. M., & Chatla, S. (n.d.). The Elephant in the Room:
#   Evaluating the Predictive Performance of Partial Least Squares (PLS) Path Models (2015).
#   SSRN Electronic Journal SSRN Journal.
#
#  And was addapted from the script:
#   https://github.com/ISS-Analytics/pls-predict/blob/master/lib/PLSpredict.R
#   in order to work directly with an plspm object.
#
# @author: Javier Lopatin
#
##########################################################################################

plspmPredict <- function(pls, testData){

    if (class(pls) != "plspm")
      stop("\n'plspmPredict()' requires a 'plspm' object")
    # test availibility of dataset (either Y or pls$data)
    plspm::test_dataset(Y, pls$data, pls$model$gens$obs)

    # =======================================================
    # inputs setting
    # =======================================================

    ltVariables = pls$model$gen$lvs_names
    mmVariables = pls$model$gen$mvs_names
    path_coef <- pls$path_coefs

    #Extract and Normalize the measurements for the model
    normDataTrain <- scale(pls$data[, mmVariables],TRUE,TRUE)

    #Extract Mean and Standard Deviation of measurements for future prediction
    meanData <- attr(normDataTrain, "scaled:center")
    sdData <- attr(normDataTrain, "scaled:scale")

    # =======================================================
    # prepare data
    # =======================================================

    #Create a matrix of outer_weights
    outer_weights <- matrix(data=0,
                            nrow=length(mmVariables),
                            ncol=length(ltVariables),
                            dimnames = list(mmVariables,ltVariables))

    # get relationship matrix
    mmMatrix = matrix(nrow = length(mmVariables),
                      ncol = 2, byrow =TRUE,
                      dimnames = list(1:length(mmVariables), c("latent","measurement")))
    mmMatrix[,'latent'] = as.character(pls$outer_model[,2])
    mmMatrix[,'measurement'] = as.character(pls$outer_model[,1])

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

    #Create a matrix of outer_loadins
    outer_loadings <- matrix(data=0,
                            nrow=length(mmVariables),
                            ncol=length(ltVariables),
                            dimnames = list(mmVariables,ltVariables))

    for (i in 1:length(ltVariables))  {
      mesVar = mmMatrix[mmMatrix[,"latent"]==ltVariables[i], "measurement"]
      idx = which(pls$outer_model$name %in% mesVar)
      outer_loadings[mesVar , ltVariables[i]] = pls$outer_model$loading[idx]
    }


    #Identify Exogenous and Endogenous Variables
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

    #Extract Measurements needed for Predictions
    normData <- testData[, pMeasurements]

    # Normalize data
    for (i in pMeasurements)
    {
      normData[,i] <-(normData[,i] - meanData[i])/sdData[i]
    }

    # Convert dataset to matrix
    normData<-data.matrix(normData)

    # Add empty columns to normData for the estimated measurements
    for (i in 1:length(eMeasurements))
    {
      normData = cbind(normData, seq(0,0,length.out =nrow(normData)))
      colnames(normData)[length(colnames(normData))]=eMeasurements[i]
    }

    # Estimate Factor Scores from Outer Path
    fscores <- normData %*% outer_weights

    # Estimate Factor Scores from Inner Path and complete Matrix
    fscores <- fscores + fscores %*% t(path_coef)

    # Predict Measurements with loadings
    predictedMeasurements <- fscores %*% t(outer_loadings)

    # Denormalize data
    for (i in mmVariables)
    {
      predictedMeasurements[,i]<-(predictedMeasurements[,i] * sdData[i])+meanData[i]
    }

    #Calculating the residuals
    residuals <- testData[,resMeasurements] - predictedMeasurements[,resMeasurements]

    #Prepare return Object
    predictResults <- list(testData = testData[,resMeasurements],
                           predictedMeasurements = predictedMeasurements[,resMeasurements],
                           residuals = residuals,
                           compositeScores = fscores)

    class(predictResults) <- "predictResults"
    return(predictResults)

}
