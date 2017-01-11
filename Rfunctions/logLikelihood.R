# Carbon balance model 
# Developed by Kashif Mahmud and Belinda Medlyn (November 2016)
# k.mahmud@westernsydney.edu.au

# This script calcualtes LogLikelihood to find the most accurate model
##############################

# Calcualte LogLi to find the most accurate model
logLikelihood <- function (data,output) {
  logLi <- matrix(0, nrow=nrow(data), ncol = 1) # Initialising the logLi
  for (i in 1:nrow(data)) {
    if (!is.na(data$Mleaf[i])) {
      logLi[i] = - 0.5*((output$Mleaf[i] - data$Mleaf[i])/data$Mleaf_SD[i])^2 - log(data$Mleaf_SD[i]) - log(2*pi)^0.5
    }
    if (!is.na(data$Mstem[i])) {
      logLi[i] = logLi[i] - 0.5*((output$Mstem[i] - data$Mstem[i])/data$Mstem_SD[i])^2 - log(data$Mstem_SD[i]) - log(2*pi)^0.5
    }
    if (!is.na(data$Mroot[i])) {
      logLi[i] = logLi[i] - 0.5*((output$Mroot[i] - data$Mroot[i])/data$Mroot_SD[i])^2 - log(data$Mroot_SD[i]) - log(2*pi)^0.5
    }
    if (!is.na(data$Sleaf[i])) {
      logLi[i] = logLi[i] - 0.5*((output$Sleaf[i] - data$Sleaf[i])/data$Sleaf_SD[i])^2 - log(data$Sleaf_SD[i]) - log(2*pi)^0.5
    }
  }
  return(sum(logLi))
}
