# plspmPred.R
This function predicts PLS-PM latent and measurement variables from a 'plspm' object ('plspm' R-package)

This is based on the publication:
   Shmueli, G., Ray, S., Estrada, J. M., & Chatla, S. (n.d.). The Elephant in the Room:
   Evaluating the Predictive Performance of Partial Least Squares (PLS) Path Models (2015).
   SSRN Electronic Journal SSRN Journal.

The script was addapted from the script:
   https://github.com/ISS-Analytics/pls-predict/blob/master/lib/PLSpredict.R
 
The adaptation was done in order to work directly with an plspm object.

Model inputs:
  - a plspm object
  - a test dataset
  
The plspmPred.R function returns a list with:
  - the test dataset
  - the prediction of all measurement variables
  - the residuals of all measurement variables
  - the predicted Latent Variable scores

@author: Javier Lopatin
