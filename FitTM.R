### Load the necessary packages. 
library('BradleyTerry2')  # Need to fit a B-T model.
library('eba')  # Need to this library to implement the pcX function (i.e.
# fit the T-M and B-T models). 
library('xtable') # Need for the xtable function. 
library('weights') # Need for the rd function. 
library('Matrix') # Need for the rankMatrix function.

FitTM <- function(data, covariates = NULL, verbose = FALSE){
  # Implement the Thurstone-Mosteller (T-M) model.
  #
  # Args: 
  #  data: Paired-comparison design matrix with pre-set row and column labels. 
  #  verbose: If TRUE, prints a summary of the model along with the
  #           95% confidence intervals for each model coefficient;
  #           if not, then nothing is printed.
  #  covariates: A matrix or data frame of additional explanatory variables
  #              (one for each column). Names of these variables is given
  #              by their corresponding column names. Defaulted to NULL.
  # Returns: 
  # The fitted T-M model fitted using the glm function. 
  #
  # Notes:
  # The pcX function creates a paired-comparison design matrix.
  # Row stimuli are chosen over column stimuli in the paired-comparison
  # design matrix.
  
  if (is.null(covariates) == TRUE){
    # Fit the T-M model using glm and compute 95% confidence intervals
    # WITHOUT additional explanatory variables. 
    m <- dim(data)[1]
    y1 <- t(data)[lower.tri(data)]
    y0 <- data[lower.tri(data)]
    TM.glm <- glm(cbind(y1, y0) ~ 0 + pcX(m, omitRef = TRUE),
                  binomial(probit))
    conf.int <- confint(TM.glm)
  } else {
    # Fit the T-M model using glm and compute 95% confidence intervals
    # WITH additional explanatory variables. 
    m <- dim(data)[1]
    y1 <- t(data)[lower.tri(data)]
    y0 <- data[lower.tri(data)]
    TM.glm <- glm(cbind(y1, y0) ~ 0 + pcX(m, omitRef = TRUE) + covariates,
                  binomial(probit))
    conf.int <- confint(TM.glm)
  }
  
  # Label the coefficients in the T-M model. 
  labels <- rownames(data)
  names(TM.glm$coefficients) <- c(labels[-1], colnames(covariates))
  
  # If verbose = TRUE, print out the results of the model along with the CIs.
  if (verbose == TRUE){
    cat('T-M Model Summary:', "\n")
    print(summary(TM.glm))
    cat("\n")
    cat('T-M 95% Confidence Intervals:', "\n")
    print(conf.int)
    cat("\n")
  }
  
  return(TM.glm)
}
