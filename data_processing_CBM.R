# Carbon balance model 
# Developed by Kashif Mahmud and Belinda Medlyn (November 2016)
# k.mahmud@westernsydney.edu.au

# This script process the raw data to carry out Bayesian calibration for 5 variables (allocation fractions: "k","Y",af","as","sf") on 
# various temporal scales to estimate Carbon pools (Cstorage,Cleaf,Cstem,Croot)
##############################

# Raw data processing for different pot volumes (e.g. 5L, 10L, ....., 1000L)
GPP.data = subset(GPP.data.raw,(volume %in% vol[v])) # Consider one pot volume at a time to run MCMC on CBM
names(GPP.data)[3] = "GPP"
Rd.data = subset(Rd.data.raw,(volume %in% vol[v]))
Sleaf.data = tnc.data = subset(tnc.data.raw,(volume %in% vol[v]))
Mleaf.data = subset(Mleaf.data.raw,(volume %in% vol[v]))
Mstem.data = subset(Mstem.data.raw,(volume %in% vol[v]))
Mroot.data = subset(Mroot.data.raw,(volume %in% vol[v]))

################ Try with inverse SE data
# Sleaf.data$tnc_SE[1:nrow(Sleaf.data)] <- Sleaf.data$tnc_SE[nrow(Sleaf.data):1]
# Mleaf.data$leafmass_SE[1:nrow(Mleaf.data)] <- Mleaf.data$leafmass_SE[nrow(Mleaf.data):1]
# Mstem.data$stemmass_SE[1:nrow(Mstem.data)] <- Mstem.data$stemmass_SE[nrow(Mstem.data):1]
# Mroot.data$rootmass_SE[1:nrow(Mroot.data)] <- Mroot.data$rootmass_SE[nrow(Mroot.data):1]

################ Try with 1/10 reduced SE data
# Sleaf.data$tnc_SE <- Sleaf.data$tnc_SE / 10
# Mleaf.data$leafmass_SE <- Mleaf.data$leafmass_SE / 10
# Mstem.data$stemmass_SE <- Mstem.data$stemmass_SE / 10
# Mroot.data$rootmass_SE <- Mroot.data$rootmass_SE / 10

# Merge all GPP, Rd, Cleaf, Cstem, Croot data
data = merge(GPP.data,Rd.data, all = TRUE)
data = merge(data,Sleaf.data, all = TRUE)
data = merge(data,Mleaf.data, all = TRUE)
data = merge(data,Mstem.data, all = TRUE)
data = merge(data,Mroot.data, all = TRUE)
names(data)[4:ncol(data)] = c("Rd","Sleaf","Sleaf_SD","Mleaf","Mleaf_SD","Mstem","Mstem_SD","Mroot","Mroot_SD")
data[ , c(7:ncol(data))] = data[ , c(7:ncol(data))] * 0.65 # Unit conversion: gDM to gC

# Reducing the measurements uncertainty (SDs) to fit the data perfectly
# data[,c("Sleaf_SD","Mleaf_SD","Mstem_SD","Mroot_SD")] = data[,c("Sleaf_SD","Mleaf_SD","Mstem_SD","Mroot_SD")] / 10000
