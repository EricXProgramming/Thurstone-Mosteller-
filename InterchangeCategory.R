### Load the necessary packages. 
library('BradleyTerry2')  # Need to fit a B-T model.
library('eba')  # Need to this library to implement the pcX function (i.e.
# fit the T-M and B-T models). 
library('xtable') # Need for the xtable function. 
library('weights') # Need for the rd function. 
library('Matrix') # Need for the rankMatrix function.

### Load functions that will be used. 
InterchangeCategory <- function(data, a, b){
  # Interchange the order of two categories of a paired-comparison 
  # design matrix with row and column labels. 
  #
  # Args: 
  #  data: Paired-comparison design matrix with row and column labels.
  #  a: Numeric index (i.e. row or column number) of first 
  #     category to be interchanged.
  #  b: Numeric index of second category to be interchanged. 
  #
  # Returns: 
  # The interchanged paired-comparison design matrix. 
  
  # Initialize some values.
  m <- dim(data)[1]
  n <- dim(data)[2]
  labels <- rownames(data)
  
  # Interchange the labels and create a new paired-comparison design matrix.    
  labels[c(a, b)] <- labels[c(b, a)] 
  new.data <- matrix(NA, nrow = m, ncol = n)
  rownames(new.data) <- labels
  colnames(new.data) <- labels
  
  # Interchange the elements.
  for(i in 1:m){
    new.index1 <- i
    if (i == a){
      new.index1 <- b
    } else if (i == b){
      new.index1 <- a
    } else {
      new.index1 <- i 
    }
    for(j in 1:n){
      if (j == a){
        new.index2 <- b
      } else if (j == b){
        new.index2 <- a
      } else {
        new.index2 <- j 
      }
      new.data[i, j] <- data[new.index1, new.index2]
    }
  }
  return(new.data)
}
