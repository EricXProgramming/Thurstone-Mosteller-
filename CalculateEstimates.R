### Load the necessary packages. 
library('BradleyTerry2')  # Need to fit a B-T model.
library('eba')  # Need to this library to implement the pcX function (i.e.
# fit the T-M and B-T models). 
library('xtable') # Need for the xtable function. 
library('weights') # Need for the rd function. 
library('Matrix') # Need for the rankMatrix function.

CalculateEstimates <- function(data, covariates = NULL){
  # Calculate three tables that consists of the estimated number of times that
  # category a is preferred over category b, estimated probability that category a 
  # is preferred over category b, and the individual pairwise deviance contributions.
  # The null deviance and model deviance is also computed. 
  
  # Args: 
  #  data: Paired-comparison design matrix with row and column labels.
  #  covariates: A matrix or data frame of additional explanatory variables
  #              (one for each column). Names of these variables is given
  #              by their corresponding column names. Defaulted to NULL.
  #
  # Returns:
  # A list containing matrices consisting of the estimated number of times that category 
  # a is preferred over category b, estimated probability that category a 
  # is preferred over category b, and the individual pairwise deviance contributions.
  # The list also contains the null deviance and model deviance. 
  
  # Initialize some objects used throughout.
  m <- dim(data)[1]
  n <- dim(data)[2]
  labels <- rownames(data)
  num.comparisons <- sum(upper.tri(data)) 
  fit.estimates <- matrix(0, nrow = m, ncol = n)  # Estimates of comparisons.
  rownames(fit.estimates) <- labels
  colnames(fit.estimates) <- labels
  true.probabilities <- matrix(0, nrow = m, ncol = n)  # True probabilities. 
  rownames(true.probabilities) <- labels
  colnames(true.probabilities) <- labels
  fit.probabilities <- matrix(0, nrow = m, ncol = n)  # Fitted probabilities. 
  rownames(fit.probabilities) <- labels
  colnames(fit.probabilities) <- labels
  deviance.contributions <- matrix(0, nrow = m, ncol = n)  # Model deviance contributions. 
  rownames(deviance.contributions) <- labels
  colnames(deviance.contributions) <- labels
  
  # Fit the T-M model using the first category as the 
  # reference category (estimate equal to 0).
  TM.glm <- FitTM(data, covariates = covariates)
  fitted.values <- fitted(TM.glm)
  
  # Compute and store all of the pairwise estimates and probablities.
  index <- 1
  for (i in 1:(m - 1)){  # Go along the rows
    for (j in (i + 1):m){  # Go along the columns. 
      pij <- fitted.values[index] 
      pji <- 1 - pij
      nij <- data[i, j]
      nji <- data[j, i]
      fit.estimates[i, j] <- pij * (nij + nji)
      fit.estimates[j, i] <- pji * (nij + nji)
      
      true.probabilities[i, j] <- nij / (nij + nji)
      true.probabilities[j, i] <- nji / (nij + nji)
      fit.probabilities[i, j] <- pij
      fit.probabilities[j, i] <- pji
      index <- index + 1
      
      # Compute the model deviance contributions. 
      devianceij1 <- nij * log(fit.probabilities[i, j])
      devianceij1[is.na(devianceij1)] <- 0
      devianceij2 <- nij * log(true.probabilities[i, j])
      devianceij2[is.na(devianceij2)] <- 0
      devianceji1 <- nji * log(fit.probabilities[j, i])
      devianceji1[is.na(devianceji1)] <- 0
      devianceji2 <- nji * log(true.probabilities[j, i])
      devianceji2[is.na(devianceji2)] <- 0
      deviance.contributions[i, j] <- -2 * (devianceij1 - devianceij2) -2 * (devianceji1 - devianceji2)
      deviance.contributions[j, i] <- 0
    }
  }
  fit.probabilities1 <- fit.probabilities[upper.tri(fit.probabilities)]
  fit.probabilities2 <- fit.probabilities[lower.tri(fit.probabilities)]
  true.probabilities1 <- true.probabilities[upper.tri(true.probabilities)]
  true.probabilities2 <- true.probabilities[lower.tri(true.probabilities)]
  true.counts1 <- data[upper.tri(data)]
  true.counts2 <- data[lower.tri(data)]
  
  # Compute the log-likelihood of the saturated model. 
  saturatedpart1 <- true.counts1 * log(true.probabilities1)
  saturatedpart1 <- saturatedpart1[!is.na(saturatedpart1)]
  saturatedpart2 <- true.counts2 * log(true.probabilities2)
  saturatedpart2 <- saturatedpart2[!is.na(saturatedpart2)]
  saturated.LL <- sum(saturatedpart1) + sum(saturatedpart2)
  
  # Compute the null deviance.
  null.deviance <- -2 * (sum(data * log(1 / 2)) - saturated.LL)
  
  # Compute the residual deviance. 
  residualpart1 <- true.counts1 * log(fit.probabilities1)
  residualpart1 <- residualpart1[!is.na(residualpart1)]
  residualpart2 <- true.counts2 * log(fit.probabilities2)
  residualpart2 <- residualpart2[!is.na(residualpart2)]
  residual.LL <- sum(residualpart1) + sum(residualpart2)
  residual.deviance <- -2 * (residual.LL - saturated.LL)
  
  # Round your answers to 2 decimal places before returning the results. 
  fit.estimates <- round(fit.estimates, 2)
  deviance.contributions <- round(deviance.contributions, 2)
  null.deviance <- round(null.deviance, 2)
  residual.deviance <- round(residual.deviance, 2)
  fit.probabilities <- round(fit.probabilities, 2)
  
  results <- list(fit.estimates, fit.probabilities, deviance.contributions, null.deviance, residual.deviance)
  names(results) <- c("Fitted Count Estimates", "Fitted Probability Estimates",
                      "Model Deviance Contributions", "Null Deviance", "Model Deviance")
  return(results)
}