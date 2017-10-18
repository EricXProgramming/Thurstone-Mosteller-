### Load the necessary packages. 
library('BradleyTerry2')  # Need to fit a B-T model.
library('eba')  # Need to this library to implement the pcX function (i.e.
# fit the T-M and B-T models). 
library('xtable') # Need for the xtable function. 
library('weights') # Need for the rd function. 
library('Matrix') # Need for the rankMatrix function.

GenerateAllComparisons <- function(data, covariates = NULL,
                                   caption = NA, printSE = FALSE){
  # Print out the main entries of a LaTeX table containing all of 
  # the paired comparisons for a Thursthone-Mosteller (T-M) model. Point
  # estimates and confidence intervals will be provided. 
  
  # Args: 
  #  data: Paired-comparison design matrix with row and column labels.
  #  caption: A string to be used for the caption of the table, defaulted
  #           to NA.
  #  printSE: If TRUE, display the standard errors (SEs) instead of the 95% CIs.
  #           Defaulted to FALSE so CIs are displayed instead.
  #  covariates: A matrix or data frame of additional explanatory variables
  #              (one for each column). Names of these variables is given
  #              by their corresponding column names. Defaulted to NULL.
  #
  # Returns:
  #  LaTeX code for a table containing the point estimates and confidence
  #  intervals or standard errors for all of the paired comparisons for 
  #  a T-M model is printed. The result is not returned as an object.
  
  # Initialize some objects used throughout.
  m <- dim(data)[1]
  n <- dim(data)[2]
  labels <- rownames(data)
  num.comparisons <- sum(upper.tri(data)) 
  estimates <- matrix(NA, nrow = m, ncol = n)  # Estimates of comparisons.
  rownames(estimates) <- labels
  colnames(estimates) <- labels
  lower.intervals <- matrix(NA, nrow = m, ncol = n)  # Lower end of the CI.
  rownames(lower.intervals) <- labels
  colnames(lower.intervals) <- labels
  upper.intervals <- matrix(NA, nrow = m, ncol = n)  # Upper end of the CI.
  rownames(upper.intervals) <- labels
  colnames(upper.intervals) <- labels
  standard.errors <- matrix(NA, nrow = m, ncol = n)  # SEs of comparisons.
  rownames(upper.intervals) <- labels
  colnames(upper.intervals) <- labels
  
  # Fit the T-M model using the first category as the
  # reference category (estimate equal to 0).
  TM.glm <- FitTM(data, covariates = covariates)
  coefficients <- coef(TM.glm)
  coefficients <- c(0, coefficients) # Add the reference category estimate. 
  covariances <- vcov(TM.glm) # No reference category.
  
  # Compute and store all of the pairwise estimates and confidence intervals.
  for (j in 1:(m - 1)){  # Go along the columns. 
    for (i in (j + 1):m){  # Go along the rows. 
      estimates[i, j] <- coefficients[i] - coefficients[j]
      if (j == 1){
        standard.errors[i, j] <- summary(TM.glm)$coef[, 2][i - 1]
      } else {
        standard.errors[i, j] <- sqrt(covariances[i - 1, i - 1] 
                                      + covariances[j - 1, j - 1] 
                                      - 2 * covariances[i - 1, j - 1])
      }
      lower.intervals[i, j] <- estimates[i, j]                                                                  - qnorm(0.975) * standard.errors[i, j]
      upper.intervals[i, j] <- estimates[i, j]                                                                  + qnorm(0.975) * standard.errors[i, j]     
    }
  }
  
  # Format and print the results into a LaTeX table.
  for (j in 1:(m - 1)){  # Go along the columns.
    if (j == 1){
      cat(paste("\\", "begin{table}[!h]", sep = ""))
      cat("\n")
      cat(paste("\\", "centering", sep = ""))
      cat("\n")
      cat(paste("\\", "caption{", caption, "}", sep = ""))
      cat("\n")
      cat(paste("\\", "begin{tabular}",
                "{ccc}", sep = ""))
      cat("\n")
      cat(paste("\t", "\\", "hline", sep = ""))
      cat("\n")
      if (printSE == TRUE){
        cat(paste("\t", "Paired Comparisons", "&", "Estimates",
                  "&", "SEs", "\\"))
        cat('\\') # Print out an extra backslash (only printed one before). 
      } else {
        cat(paste("\t", "Paired Comparisons", "&", "Estimates",
                  "&", "95\\% C.I.", "\\"))
        cat('\\') # Print out an extra backslash (only printed one before). 
      }
      cat("\n")
      cat(paste("\t", "\\", "hline", sep = ""))
      cat("\n")
    }
    for(i in (j + 1):m){  # Go along the rows.
      category1 <- labels[i]
      category2 <- labels[j]
      estimate <- rd(estimates[i, j], 3)
      lower.interval <- rd(lower.intervals[i, j], 3)
      upper.interval <- rd(upper.intervals[i, j], 3)
      standard.error <- rd(standard.errors[i, j], 3)
      if (printSE == TRUE){
        cat(paste("\t", category1, "vs.", category2,
                  "&", estimate, "&", standard.error, "\\"))
        cat('\\') # Print out an extra backslash (only printed one before).
      } else {
        cat(paste("\t", category1, "vs.", category2,
                  "&", estimate,
                  "&", "$(", lower.interval, ",", upper.interval, ")$",
                  "\\"))
        cat('\\') # Print out an extra backslash (only printed one before).
      }
      cat("\n")
    }
    cat(paste("\t", "\\", "hline", sep = ""))
    cat("\n")
    if(j == (m - 1)){
      cat(paste("\\", "end{tabular}", sep = ""))
      cat("\n")
      cat(paste("\\", "end{table}", sep = ""))
    }
  }
}