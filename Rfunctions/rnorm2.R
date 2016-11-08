# Function to gerenate normally distributed random numbers with mean and SD

rnorm2 <- function(n,mean,sd) { 
  mean+sd*scale(rnorm(n)) 
  }