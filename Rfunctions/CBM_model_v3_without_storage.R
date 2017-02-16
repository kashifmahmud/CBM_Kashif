# Carbon balance model 
# Developed by Kashif Mahmud and Belinda Medlyn (November 2016)
# k.mahmud@westernsydney.edu.au

# This script define the model equations to carry out Bayesian calibration for 5 variables (allocation fractions: "k","Y",af","as","sf") on 
# various temporal scales to estimate Carbon pools (Cstorage,Cleaf,Cstem,Croot)
##############################
# This version is for previous model that does not consider storage pool

# Defining the model to iteratively calculate Cstorage, Cleaf, Cstem, Croot, Sleaf, Sstem, Sroot
model <- function (GPP,Rd,no.param,Mleaf,Mstem,Mroot,Y,af,as,sf) {
  
  if (no.param == 1) {
    Y.i = Y[1]; af.i = af[1]; as.i = as[1]; sf.i = sf[1]
  }
  
  for (i in 2:length(GPP)) {
    if (no.param == 2) {
      Y.i = Y[1]+ Y[2]*i; af.i = af[1]+ af[2]*i; as.i = as[1]+ as[2]*i; sf.i = sf[1]+ sf[2]*i
      # sf.i = 0
      # if (i > length(GPP)/2) {
      # sf.i = sf[1]+ sf[2]*i
      # }
    }
    if (no.param == 3) {
      Y.i = Y[1]+ Y[2]*i + Y[3]*i*i; af.i = af[1]+ af[2]*i + af[3]*i*i; 
      as.i = as[1]+ as[2]*i + as[3]*i*i; 
      sf.i = sf[1]+ sf[2]*i + sf[3]*i*i
      # sf.i = 0
      # if (i > length(GPP)/2) {
      #   sf.i = sf[1]+ sf[2]*i + sf[3]*i*i
      # }
    }
    
    Mleaf[i] <- Mleaf[i-1] + GPP[i-1]*af.i*(1-Y.i) - sf.i*Mleaf[i-1]
    Mstem[i] <- Mstem[i-1] + GPP[i-1]*as.i*(1-Y.i)
    Mroot[i] <- Mroot[i-1] + GPP[i-1]*(1-af.i-as.i)*(1-Y.i)
  }
  output = data.frame(Mleaf,Mstem,Mroot)
  # print("successfully calculated Carbon pools")
  
  return(output)
}
