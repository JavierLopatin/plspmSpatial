# plspmPred.R
This function predicts PLS-PM latent and measurement variables from a 'plspm' object ('plspm' R-package)

This is based on the publication:
   Shmueli, G., Ray, S., Estrada, J. M., & Chatla, S. (n.d.). The Elephant in the Room:
   Evaluating the Predictive Performance of Partial Least Squares (PLS) Path Models (2015).
   SSRN Electronic Journal SSRN Journal.

The script was adapted from the script:
   https://github.com/ISS-Analytics/pls-predict/blob/master/lib/PLSpredict.R

The adaptation was done in order to work directly with an plspm object.

WARNING: For the moment, only working with <code>class(data) == data.frame</code> object for prediction. Raster classes and accuracy metrics to be added

------------------------

Usage:
-----
  plspmPred(pls, data, ...)

Arguments:
-----------
   - pls: An plspm object from the plspm package
   - dat: dataframe or Raster Stack with the model predictors

Details:
-------
The function plspmPredict estimates  extrapolation values

Values:
------
An object of class <code>plspmPredict</code> is returned. The object returns a list with:
  - **mmPredicted**
        data.frame or RasterStack of the predicted all measurement variables
  - **mmResiduals**
        data.frame or RasterStack of the residuals of all measurement variables
  - **lvPredict**
        data.frame or RasterStack of the predicted Latent Variables scores [in ordination units]
  - **lvResiduals**
        data.frame or RasterStack of the residuals of all Latent Variables scores

@author: Javier Lopatin
