# Carbon balance model parameter sensitivity analysis
# Developed by Kashif Mahmud and Belinda Medlyn (November 2016)
# k.mahmud@westernsydney.edu.au

# This script define the model equations to carry out Bayesian calibration for 5 variables (allocation fractions: "k","Y",af","as","sf") on 
# various temporal scales to estimate Carbon pools (Cstorage,Cleaf,Cstem,Croot)
##############################

# Defining the model to iteratively calculate Cstorage, Cleaf, Cstem, Croot, Sleaf, Sstem, Sroot
model_v2 <- function (Cday,Rd,no.param,Mleaf,Mstem,Mroot,leaf.data,b,intercept,k,Y,af,as,sf) {
  GPP = Cstorage = Sleaf = Sstem = Sroot = M = c()
  
  # From Duan's experiment for TNC partitioning to tree organs
  # Leaf TNC/Leaf DW =  0.1401421; Stem TNC/Stem DW =  0.0453869; Root TNC/Root DW =  0.02154037
  Sleaf[1] = Mleaf[1] / 0.65 * 0.1401421
  Sstem[1] = Mstem[1] / 0.65 * 0.0453869
  Sroot[1] = Mroot[1] / 0.65 * 0.02154037
  Cstorage[1] <- Sleaf[1] + Sstem[1] + Sroot[1] 
  
  Cleaf <- Croot <- Cstem <- LA <- c()
  Cleaf[1] <- Mleaf[1] - Sleaf[1]
  Cstem[1] <- Mstem[1] - Sstem[1]
  Croot[1] <- Mroot[1] - Sroot[1]
  
  LA[1] <- leaf.data$initial_LA
  
  if (no.param == 1) {
    k.i = k[1]; Y.i = Y[1]; af.i = af[1]; as.i = as[1]; sf.i = sf[1]
  }
  
  for (i in 2:length(Cday)) {
  #   if (no.param == 2) {
  #     k.i = k[1] + k[2]*i; Y.i = Y[1]+ Y[2]*i; af.i = af[1]+ af[2]*i; as.i = as[1]+ as[2]*i; sf.i = sf[1]+ sf[2]*i
  #     # sf.i = 0
  #     # if (i > length(GPP)/2) {
  #     # sf.i = sf[1]+ sf[2]*i
  #     # }
  #   }
  #   if (no.param == 3) {
  #     k.i = k[1] + k[2]*i + k[3]*i*i; Y.i = Y[1]+ Y[2]*i + Y[3]*i*i; af.i = af[1]+ af[2]*i + af[3]*i*i; 
  #     as.i = as[1]+ as[2]*i + as[3]*i*i; 
  #     sf.i = sf[1]+ sf[2]*i + sf[3]*i*i
  #     # sf.i = 0
  #     # if (i > length(GPP)/2) {
  #     #   sf.i = sf[1]+ sf[2]*i + sf[3]*i*i
  #     # }
  #   }
    
    M[i-1] <- b * LA[i-1] + intercept
    GPP[i-1] <- LA[i-1] * Cday[i-1] * M[i-1] # calculate total daily C gain with self shading
    
    
    Cstorage[i] <- Cstorage[i-1] + GPP[i-1] - Rd[i-1]*(Mleaf[i-1] + Mroot[i-1] + Mstem[i-1]) - k.i*Cstorage[i-1]
    Sleaf[i] <- Cstorage[i] * 0.75 # 75% of storage goes to leaf (Duan's experiment)
    Sstem[i] <- Cstorage[i] * 0.16 # 16% of storage goes to stem (Duan's experiment)
    Sroot[i] <- Cstorage[i] * 0.09 # 9% of storage goes to root (Duan's experiment)
    
    Cleaf[i] <- Cleaf[i-1] + k.i*Cstorage[i-1]*af.i*(1-Y.i) - sf.i*Cleaf[i-1]
    Cstem[i] <- Cstem[i-1] + k.i*Cstorage[i-1]*as.i*(1-Y.i)
    Croot[i] <- Croot[i-1] + k.i*Cstorage[i-1]*(1-af.i-as.i)*(1-Y.i)
    
    Mleaf[i] <- Cleaf[i] + Sleaf[i]
    Mstem[i] <- Cstem[i] + Sstem[i]
    Mroot[i] <- Croot[i] + Sroot[i]
    
    # Leaf area (t) = Leaf area (T) * Leaf count (t) / Leaf count (T); t = time, T = time of harvest
    LA[i] <- leaf.data$final_LA * Mleaf[i] / leaf.data$final_LM
  }
  output = data.frame(Cstorage,Mleaf,Mstem,Mroot,Sleaf)
  # print("successfully calculated Carbon pools")
  
  return(output)
}
