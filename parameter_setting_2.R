

# Import daily GPP, daily Rd data
# c3 = read.csv("rawdata/c3.csv") # min and max of c3 coefficients with quardatic equation fit for all paramter set

# Setting lower and upper bounds of the prior parameter pdf, and starting point of the chain
param.k <- matrix(c(0,0.45,1) , nrow=1, ncol=3, byrow=T) 
param.Y <- matrix(c(0.25,0.3,0.35) , nrow=1, ncol=3, byrow=T) 
param.af <- matrix(c(0.1,0.4,0.7) , nrow=1, ncol=3, byrow=T) 
param.as <- matrix(c(0.1,0.25,0.6) , nrow=1, ncol=3, byrow=T) 
# param.sf <- matrix(c(0,0.05,0.1) , nrow=1, ncol=3, byrow=T)

if (no.param > 1) {
  param.k <- rbind(param.k, c(-(param.k[3]-param.k[1])/nrow(data), 0, (param.k[3]-param.k[1])/nrow(data)))
  param.Y <- rbind(matrix(c(0.28,0.3,0.32) , nrow=1, ncol=3, byrow=T), c(-(param.Y[3]-param.Y[1])/2/nrow(data), 0, (param.Y[3]-param.Y[1])/2/nrow(data)))
  param.af <- rbind(param.af, c(-(param.af[3]-param.af[1])/nrow(data), 0, (param.af[3]-param.af[1])/nrow(data)))
  param.as <- rbind(param.as, c(-(param.as[3]-param.as[1])/2/nrow(data), 0, (param.as[3]-param.as[1])/2/nrow(data)))
  # param.sf <- rbind(param.sf, c(-(param.sf[3]-param.sf[1])/nrow(data), 0, (param.sf[3]-param.sf[1])/nrow(data)))
}

if (no.param > 2) {
  # param.k <- rbind(param.k, c(c3$k[1], 0, c3$k[2]))
  # param.Y <- rbind(param.Y, c(c3$Y[1], 0, c3$Y[2]))
  # param.af <- rbind(param.af, c(c3$af[1], 0, c3$af[2]))
  # param.as <- rbind(param.as, c(c3$as[1], 0, c3$as[2]))
  # param.sf <- rbind(param.sf, c(c3$sf[1], 0, c3$sf[2]))
  param.k <- rbind(param.k, c((param.k[1,1]-param.k[1,3]-param.k[2,3]*nrow(data))/(nrow(data)^2), 0, (param.k[1,3]-param.k[1,1]-param.k[2,1]*nrow(data))/(nrow(data)^2)))
  param.Y <- rbind(param.Y, c((param.Y[1,1]-param.Y[1,3]-param.Y[2,3]*nrow(data))/(2*nrow(data)^2), 0, (param.Y[1,3]-param.Y[1,1]-param.Y[2,1]*nrow(data))/(2*nrow(data)^2)))
  param.af <- rbind(param.af, c((param.af[1,1]-param.af[1,3]-param.af[2,3]*nrow(data))/(nrow(data)^2), 0, (param.af[1,3]-param.af[1,1]-param.af[2,1]*nrow(data))/(nrow(data)^2)))
  param.as <- rbind(param.as, c((param.as[1,1]-param.as[1,3]-param.as[2,3]*nrow(data))/(2*nrow(data)^2), 0, (param.as[1,3]-param.as[1,1]-param.as[2,1]*nrow(data))/(2*nrow(data)^2)))
  # param.sf <- rbind(param.sf, c((param.sf[1,1]-param.sf[1,3]-param.sf[2,3]*nrow(data))/(nrow(data)^2), 0, (param.sf[1,3]-param.sf[1,1]-param.sf[2,1]*nrow(data))/(nrow(data)^2)))
}

param = data.frame(param.k,param.Y,param.af,param.as)
names(param) <- c("k_min","k","k_max","Y_min","Y","Y_max","af_min","af","af_max","as_min","as","as_max")
pMinima <- param[ ,c("k_min","Y_min","af_min","as_min")]
pMaxima <- param[ ,c("k_max","Y_max","af_max","as_max")]
pValues <- param[ ,c("k","Y","af","as")] # Starting point of the chain
pChain <- matrix(0, nrow=chainLength, ncol = no.param*no.var+1) # Initialising the chain

