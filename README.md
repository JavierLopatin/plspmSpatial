# plspmTools

Set of functions to help in PLS-PM analysis of ecological and geoscience data.

The currently included functions are:

    - plspmPredict
    - plspm.groupsPredict
    - plspmResiduals

@author: Javier Lopatin

### **plspmPredict**:

This function predicts PLS-PM latent and measurement variables from a 'plspm' object ('plspm' R-package)

This is based on the publication:
   Shmueli, G., Ray, S., Estrada, J. M., & Chatla, S. (n.d.). The Elephant in the Room:
   Evaluating the Predictive Performance of Partial Least Squares (PLS) Path Models (2015).
   SSRN Electronic Journal SSRN Journal.

The script was adapted from the script:
   <https://github.com/ISS-Analytics/pls-predict/blob/master/lib/PLSpredict.R>

The adaptation was done in order to work directly with an plspm object.

WARNING: For the moment, only working with <code>class(data) == data.frame</code> object for prediction. Raster classes to be added


**Usage**:

  plspmPred(pls, dat, ...)

**Arguments**:

-   pls: An plspm object from the plspm package
-   dat: data.frame or Raster Stack with the model predictors

**Details**:

The function plspmPredict estimates  extrapolation values of Latent and Measurement Variables from and plspm object

**Values**:

An object of class <code>plspmPredict</code> is returned. The object returns a list with:

-   **mmData**

        Matrix or RasterStack of the input measurement variables

-   **mmPredicted**

        Matrix or RasterStack of the predicted all measurement variables

-   **mmResiduals**

        Matrix or RasterStack of the residuals of all measurement variables

-   **Scores**

        Matrix or RasterStack of the predicted Latent Variables scores [in ordination units]

-   **r_square**

        Matrix of Squared Pearson's Correlation values of all measurement variables (only with 'dat' as class 'data.frame')

-   **RMSE**

        Matrix of Root-Mean-Square-Error values of all measurement variables (only with 'dat' as class 'data.frame')

-   **nRMSE**

        Matrix of normalizedRoot-Mean-Square-Error [%] values of all measurement variables (only with 'dat' as class 'data.frame')

-   **bias**

        Matrix of bias values of all measurement variables (only with 'dat' as class 'data.frame')

### **plspm.groupsPredict**:

**Usage**:

  plspm.groupsPredict(pls, pls.groups, train.groups, dat)

  **Details**:

This function has the same functions as <code>plspmPredict</code>, but uses a <code>plspm.groups</code> object as input. Therefore, it gives a list of predicted scores, measurement variables, and residuals for the 'general', 'group1', and 'group2' models.

### **plspmRsiduals**:

**Usage**:

  plspmRsiduals(pls)

**Arguments**:

-   pls: An plspm object from the plspm package
-   dat: data.frame or Raster Stack with the model predictors

**Details**:

This function obtain residuals for all Latent and Measurement variables from an 'plspm' object

**Values**:

An object of class <code>plspmResiduals</code> is returned. The object returns a list with:

-   **inner_residuals**

        Matrix of residual values for the Latent Variables

-   **outer_residuals**

        Matrix of residual values for the Measurement Variables
