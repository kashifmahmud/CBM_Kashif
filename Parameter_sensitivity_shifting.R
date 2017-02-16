# Carbon balance model (CBM)
# Developed by Kashif Mahmud and Belinda Medlyn (November 2016)
# k.mahmud@westernsydney.edu.au

# This code carries out Parameter sensitivity analyses for 6 variables (allocation fractions: "k","Y",af","as","ar","sf") 
# to estimate Carbon pools (Cstorage,Cleaf,Cstem,Croot)

##############################
# Set working directory for saving figures
setwd("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif")

# This script cleans the workspace, loads necessary Rfunctions and packages
source("load_packages_functions_CBM.R")

# # Load the function to define the CBM equations to iteratively calculate Cstorage, Cleaf, Cstem, Croot, Sleaf, Sstem, Sroot
# source("Rfunctions/CBM_model_v2.R")

# Import daily Cday, daily Rd data
Cday.data.raw <- read.csv("rawdata/cday_120_clean_gross.csv") # Units gC d-1
Cday.data.raw$Date <- as.Date(Cday.data.raw$Date)
# Cday.data.raw = read.csv("rawdata/GPP.csv")
# GPP.data.raw = read.csv("rawdata/GPP.csv") # Units gC d-1
Rd.data.raw = read.csv("rawdata/Rd.csv") # Units g C g-1 plant d-1
tnc.data.raw = read.csv("rawdata/tnc_fortnightly_data.csv") # Units g plant-1

# Import weekly Cleaf, weekly Cstem, initial/harvest Croot data with Mean and SD
Mleaf.data.raw = read.csv("rawdata/Cleaf_weekly_data.csv") # Units gC d-1
Mstem.data.raw = read.csv("rawdata/Cstem_weekly_data.csv") # Units gC d-1
Mroot.data.raw = read.csv("rawdata/Croot_twice_data.csv") # Units gC d-1

# Import daily modelled parameters
param.summary = read.csv("output/processeddata/summary.param.csv")

################### Harvested leaf mass and leaf area for 5L pot
plot.summary = read.csv("rawdata/plot_summary.csv")
harvest_data = read.csv("rawdata/harvest aboveground mass.csv")
harvest_data <- merge(plot.summary,harvest_data,by=c("pot","plot"))
hd.idn = subset(harvest_data,volume==5) 
hd.idn[nrow(hd.idn)+1, ] = colMeans(hd.idn, na.rm = TRUE) # R8 = Average of leaf counts
hd.idn[nrow(hd.idn)+1, ] = (apply(hd.idn, 2, sd))/(nrow(hd.idn)-1)^0.5 # R9 = Standard deviation of leaf counts
hd.idn$newleaf_count <- round(hd.idn$newleaf_count)
hd.idn$leaf_count <- round(hd.idn$leaf_count)
hd.idn$leaf_mass <- hd.idn$Leafmass.bag - hd.idn$Leaf_bag
hd.idn$stem_mass <- hd.idn$stemmass.bag - hd.idn$stem_bag
keeps <- c("volume", "leaf_area", "leaf_count", "leaf_mass", "stem_mass", "Croot", "Froot")
hd = hd.idn[ , keeps, drop = FALSE]
dimnames(hd)[[1]] <- c(1:7, "mean", "SE")
hd.final <- hd[0,]
hd.final <- hd["mean", ]
hd.final$leaf_area <- hd.final$leaf_area / (100*100) # unit conversion from cm2 to m2

# Import initial seedling data
initial.data <- read.csv("rawdata/seedling_initial.csv")
initial.data[nrow(initial.data)+1, 2:ncol(initial.data)] = colMeans(initial.data[2:ncol(initial.data)], na.rm = TRUE) # R7 = Average of leaf data
initial.data[nrow(initial.data)+1, 2:ncol(initial.data)] = (apply(initial.data[2:ncol(initial.data)], 2, sd, na.rm = TRUE))/(nrow(initial.data)-1)^0.5 # R8 = Standard deviation of leaf counts
# initial.data[nrow(initial.data)+1, 2:ncol(initial.data)] = apply(initial.data[2:ncol(initial.data)], 2, sd, na.rm = TRUE) # R8 = Standard deviation of leaf counts
dimnames(initial.data)[[1]] <- c(1:(nrow(initial.data)-2), "mean", "SE")
id.final <- initial.data[0,]
id.final <- initial.data["mean", ]
id.final$leaf_numb <- round(id.final$leaf_numb)
# Leaf area (at t=1) = Leaf area (T) * Leaf count (t=1) / Leaf count (T); t = 1 (initial), T = time of harvest
id.final$leaf_area = hd.final$leaf_area * id.final$leaf_numb / hd.final$leaf_count

# Self shading factors for 5L pot (use slope intercept, from 5-free vol)
sigma <- read.csv("rawdata/M_leafarea_model.csv")
sigma = subset(sigma,volume==5)
b = sigma$b
intercept = sigma$intercept

# Raw data processing for 5L pot only (1000L)
Cday.data = subset(Cday.data.raw,(volume %in% 15)) # Consider the free seedling to test the parameter sensitivity
names(Cday.data)[3] = "Cday"
Rd.data = subset(Rd.data.raw,volume==15)
Sleaf.data = tnc.data = subset(tnc.data.raw,volume==15)
Mleaf.data = subset(Mleaf.data.raw,volume==15)
Mstem.data = subset(Mstem.data.raw,volume==15)
Mroot.data = subset(Mroot.data.raw,volume==15)


# Merge all GPP, Rd, Cleaf, Cstem, Croot data
data = merge(Cday.data,Rd.data, all = TRUE)
data = merge(data,Sleaf.data, all = TRUE)
data = merge(data,Mleaf.data, all = TRUE)
data = merge(data,Mstem.data, all = TRUE)
data = merge(data,Mroot.data, all = TRUE)
names(data)[4:ncol(data)] = c("Rd","Sleaf","Sleaf_SD","Mleaf","Mleaf_SD","Mstem","Mstem_SD","Mroot","Mroot_SD")
data[ , c(7:ncol(data))] = data[ , c(7:ncol(data))] * 0.65 # Unit conversion: gDM to gC

# Form data frame with necessary leaf data
leaf.data = data.frame(id.final$leaf_numb, id.final$leaf_area, hd.final$leaf_count, hd.final$leaf_area, Mleaf.data$leafmass[nrow(Mleaf.data)])
names(leaf.data)[1:ncol(leaf.data)] <- c("initial_LC", "initial_LA", "final_LC", "final_LA", "final_LM")

##################------------------------------
# Take the parameter values for free seedling (1000L)
q=1 # Case 1
param = subset(param.summary,(volume.group %in% 2))
keeps = c("Date", "variable", "Parameter")
param = param[ , keeps, drop = FALSE]
param.casted = dcast( param , Date ~ variable )

# Calculating model outputs
Mleaf = Mstem = Mroot = c()
Mleaf[1] <- data$Mleaf[1]
Mstem[1] <- data$Mstem[1]
Mroot[1] <- data$Mroot[1]

Cday=data$Cday; Rd=data$Rd;
k=param.casted$k; Y=param.casted$Y; af=param.casted$af; as=param.casted$as; ar=param.casted$ar; sf=param.casted$sf

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

for (i in 2:length(Cday)) {
  M[i-1] <- b * LA[i-1] + intercept
  GPP[i-1] <- LA[i-1] * Cday[i-1] * M[i-1] # calculate total daily C gain with self shading
  
  Cstorage[i] <- Cstorage[i-1] + GPP[i-1] - Rd[i-1]*(Mleaf[i-1] + Mroot[i-1] + Mstem[i-1]) - k*Cstorage[i-1]
  Sleaf[i] <- Cstorage[i] * 0.75 # 75% of storage goes to leaf (Duan's experiment)
  Sstem[i] <- Cstorage[i] * 0.16 # 16% of storage goes to stem (Duan's experiment)
  Sroot[i] <- Cstorage[i] * 0.09 # 9% of storage goes to root (Duan's experiment)
  
  Cleaf[i] <- Cleaf[i-1] + k*Cstorage[i-1]*af*(1-Y) - sf*Cleaf[i-1]
  Cstem[i] <- Cstem[i-1] + k*Cstorage[i-1]*as*(1-Y)
  Croot[i] <- Croot[i-1] + k*Cstorage[i-1]*(1-af-as)*(1-Y)
  
  Mleaf[i] <- Cleaf[i] + Sleaf[i]
  Mstem[i] <- Cstem[i] + Sstem[i]
  Mroot[i] <- Croot[i] + Sroot[i]
  
  # Leaf area (t) = Leaf area (T) * Leaf count (t) / Leaf count (T); t = time, T = time of harvest
  LA[i] <- leaf.data$final_LA * Cleaf[i] / leaf.data$final_LM
}
output.final = data.frame(Cstorage,Mleaf,Mstem,Mroot,Sleaf)

# Plant Carbon pools for various parameter sensitivity
output.final$Date = data$Date
names(output.final) = c("Cstorage","Mleaf","Mstem","Mroot","Sleaf","Date")
melted.output = melt(output.final[,c("Mleaf","Mstem","Mroot","Cstorage","Sleaf","Date")], id.vars="Date")
melted.Cstorage = output.final[,c("Cstorage","Date")]
melted.output$Date = as.Date(melted.output$Date)
melted.Cstorage$Date = as.Date(melted.Cstorage$Date)
melted.output$Case = as.factor(q)
melted.Cstorage$Case = as.factor(q)

# Storing the summary of data, outputs, Cstorage, parameters
if (q == 1) {
  summary.output = melted.output
  summary.Cstorage = melted.Cstorage
}
if (q > 1) {
  summary.output = rbind(summary.output,melted.output)
  summary.Cstorage = rbind(summary.Cstorage,melted.Cstorage)
}

# p1 = ggplot() +
#   geom_line(data = summary.output, aes(x = Date, y = value, group = variable, colour=variable)) +
#   geom_point(size=2) +
#   ylab("C pools") +
#   ggtitle(paste("C pools")) +
#   # labs(colour="C pools", linetype=as.character(var[p]), shape=as.character(var[p])) +
#   theme_bw() +
#   theme(plot.title = element_text(size = 12, face = "bold")) +
#   theme(legend.title = element_text(colour="chocolate", size=12, face="bold")) +
#   theme(axis.title.x = element_text(size = 12, vjust=-.2)) +
#   theme(axis.title.y = element_text(size = 12, vjust=0.3))
# p1
# ggsave(p1,filename=paste("output/figures/Carbon_pools_case_",q,".png",sep=""))


############----------------------------------------
# Take the parameters af, as, ar for 5L pot
q=2 # Case 2
param.free = subset(param.summary,(volume.group %in% 2 & variable %in% c("k","Y","sf")))
param.5L = subset(param.summary,(volume.group %in% 1 & variable %in% c("af","as","ar")))
param = rbind(param.free, param.5L)
keeps = c("Date", "variable", "Parameter")
param = param[ , keeps, drop = FALSE]
param.casted = dcast( param , Date ~ variable )

# Calculating model outputs
Mleaf = Mstem = Mroot = c()
Mleaf[1] <- data$Mleaf[1]
Mstem[1] <- data$Mstem[1]
Mroot[1] <- data$Mroot[1]

Cday=data$Cday; Rd=data$Rd;
k=param.casted$k; Y=param.casted$Y; af=param.casted$af; as=param.casted$as; ar=param.casted$ar; sf=param.casted$sf

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

for (i in 2:length(Cday)) {
  M[i-1] <- b * LA[i-1] + intercept
  GPP[i-1] <- LA[i-1] * Cday[i-1] * M[i-1] # calculate total daily C gain with self shading
  
  Cstorage[i] <- Cstorage[i-1] + GPP[i-1] - Rd[i-1]*(Mleaf[i-1] + Mroot[i-1] + Mstem[i-1]) - k*Cstorage[i-1]
  Sleaf[i] <- Cstorage[i] * 0.75 # 75% of storage goes to leaf (Duan's experiment)
  Sstem[i] <- Cstorage[i] * 0.16 # 16% of storage goes to stem (Duan's experiment)
  Sroot[i] <- Cstorage[i] * 0.09 # 9% of storage goes to root (Duan's experiment)
  
  Cleaf[i] <- Cleaf[i-1] + k*Cstorage[i-1]*af*(1-Y) - sf*Cleaf[i-1]
  Cstem[i] <- Cstem[i-1] + k*Cstorage[i-1]*as*(1-Y)
  Croot[i] <- Croot[i-1] + k*Cstorage[i-1]*(1-af-as)*(1-Y)
  
  Mleaf[i] <- Cleaf[i] + Sleaf[i]
  Mstem[i] <- Cstem[i] + Sstem[i]
  Mroot[i] <- Croot[i] + Sroot[i]
  
  # Leaf area (t) = Leaf area (T) * Leaf count (t) / Leaf count (T); t = time, T = time of harvest
  LA[i] <- leaf.data$final_LA * Cleaf[i] / leaf.data$final_LM
}
output.final = data.frame(Cstorage,Mleaf,Mstem,Mroot,Sleaf)

# Plant Carbon pools for various parameter sensitivity
output.final$Date = data$Date
names(output.final) = c("Cstorage","Mleaf","Mstem","Mroot","Sleaf","Date")
melted.output = melt(output.final[,c("Mleaf","Mstem","Mroot","Cstorage","Sleaf","Date")], id.vars="Date")
melted.Cstorage = output.final[,c("Cstorage","Date")]
melted.output$Date = as.Date(melted.output$Date)
melted.Cstorage$Date = as.Date(melted.Cstorage$Date)
melted.output$Case = as.factor(q)
melted.Cstorage$Case = as.factor(q)

# Storing the summary of data, outputs, Cstorage, parameters
if (q == 1) {
  summary.output = melted.output
  summary.Cstorage = melted.Cstorage
}
if (q > 1) {
  summary.output = rbind(summary.output,melted.output)
  summary.Cstorage = rbind(summary.Cstorage,melted.Cstorage)
}

############----------------------------------------
# Take the parameters Y, af, as, ar for 5L pot
q=3 # Case 3
param.free = subset(param.summary,(volume.group %in% 2 & variable %in% c("k","sf")))
param.5L = subset(param.summary,(volume.group %in% 1 & variable %in% c("Y","af","as","ar")))
param = rbind(param.free, param.5L)
keeps = c("Date", "variable", "Parameter")
param = param[ , keeps, drop = FALSE]
param.casted = dcast( param , Date ~ variable )

# Calculating model outputs
Mleaf = Mstem = Mroot = c()
Mleaf[1] <- data$Mleaf[1]
Mstem[1] <- data$Mstem[1]
Mroot[1] <- data$Mroot[1]

Cday=data$Cday; Rd=data$Rd;
k=param.casted$k; Y=param.casted$Y; af=param.casted$af; as=param.casted$as; ar=param.casted$ar; sf=param.casted$sf

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

for (i in 2:length(Cday)) {
  M[i-1] <- b * LA[i-1] + intercept
  GPP[i-1] <- LA[i-1] * Cday[i-1] * M[i-1] # calculate total daily C gain with self shading
  
  Cstorage[i] <- Cstorage[i-1] + GPP[i-1] - Rd[i-1]*(Mleaf[i-1] + Mroot[i-1] + Mstem[i-1]) - k*Cstorage[i-1]
  Sleaf[i] <- Cstorage[i] * 0.75 # 75% of storage goes to leaf (Duan's experiment)
  Sstem[i] <- Cstorage[i] * 0.16 # 16% of storage goes to stem (Duan's experiment)
  Sroot[i] <- Cstorage[i] * 0.09 # 9% of storage goes to root (Duan's experiment)
  
  Cleaf[i] <- Cleaf[i-1] + k*Cstorage[i-1]*af*(1-Y) - sf*Cleaf[i-1]
  Cstem[i] <- Cstem[i-1] + k*Cstorage[i-1]*as*(1-Y)
  Croot[i] <- Croot[i-1] + k*Cstorage[i-1]*(1-af-as)*(1-Y)
  
  Mleaf[i] <- Cleaf[i] + Sleaf[i]
  Mstem[i] <- Cstem[i] + Sstem[i]
  Mroot[i] <- Croot[i] + Sroot[i]
  
  # Leaf area (t) = Leaf area (T) * Leaf count (t) / Leaf count (T); t = time, T = time of harvest
  LA[i] <- leaf.data$final_LA * Cleaf[i] / leaf.data$final_LM
}
output.final = data.frame(Cstorage,Mleaf,Mstem,Mroot,Sleaf)

# Plant Carbon pools for various parameter sensitivity
output.final$Date = data$Date
names(output.final) = c("Cstorage","Mleaf","Mstem","Mroot","Sleaf","Date")
melted.output = melt(output.final[,c("Mleaf","Mstem","Mroot","Cstorage","Sleaf","Date")], id.vars="Date")
melted.Cstorage = output.final[,c("Cstorage","Date")]
melted.output$Date = as.Date(melted.output$Date)
melted.Cstorage$Date = as.Date(melted.Cstorage$Date)
melted.output$Case = as.factor(q)
melted.Cstorage$Case = as.factor(q)

# Storing the summary of data, outputs, Cstorage, parameters
if (q == 1) {
  summary.output = melted.output
  summary.Cstorage = melted.Cstorage
}
if (q > 1) {
  summary.output = rbind(summary.output,melted.output)
  summary.Cstorage = rbind(summary.Cstorage,melted.Cstorage)
}

############----------------------------------------
# Take the parameters Y, af, as, ar, sf for 5L pot
q=4 # Case 4
param.free = subset(param.summary,(volume.group %in% 2 & variable %in% c("k")))
param.5L = subset(param.summary,(volume.group %in% 1 & variable %in% c("Y","af","as","ar","sf")))
param = rbind(param.free, param.5L)
keeps = c("Date", "variable", "Parameter")
param = param[ , keeps, drop = FALSE]
param.casted = dcast( param , Date ~ variable )

# Calculating model outputs
Mleaf = Mstem = Mroot = c()
Mleaf[1] <- data$Mleaf[1]
Mstem[1] <- data$Mstem[1]
Mroot[1] <- data$Mroot[1]

Cday=data$Cday; Rd=data$Rd;
k=param.casted$k; Y=param.casted$Y; af=param.casted$af; as=param.casted$as; ar=param.casted$ar; sf=param.casted$sf

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

for (i in 2:length(Cday)) {
  M[i-1] <- b * LA[i-1] + intercept
  GPP[i-1] <- LA[i-1] * Cday[i-1] * M[i-1] # calculate total daily C gain with self shading
  
  Cstorage[i] <- Cstorage[i-1] + GPP[i-1] - Rd[i-1]*(Mleaf[i-1] + Mroot[i-1] + Mstem[i-1]) - k*Cstorage[i-1]
  Sleaf[i] <- Cstorage[i] * 0.75 # 75% of storage goes to leaf (Duan's experiment)
  Sstem[i] <- Cstorage[i] * 0.16 # 16% of storage goes to stem (Duan's experiment)
  Sroot[i] <- Cstorage[i] * 0.09 # 9% of storage goes to root (Duan's experiment)
  
  Cleaf[i] <- Cleaf[i-1] + k*Cstorage[i-1]*af*(1-Y) - sf*Cleaf[i-1]
  Cstem[i] <- Cstem[i-1] + k*Cstorage[i-1]*as*(1-Y)
  Croot[i] <- Croot[i-1] + k*Cstorage[i-1]*(1-af-as)*(1-Y)
  
  Mleaf[i] <- Cleaf[i] + Sleaf[i]
  Mstem[i] <- Cstem[i] + Sstem[i]
  Mroot[i] <- Croot[i] + Sroot[i]
  
  # Leaf area (t) = Leaf area (T) * Leaf count (t) / Leaf count (T); t = time, T = time of harvest
  LA[i] <- leaf.data$final_LA * Cleaf[i] / leaf.data$final_LM
}
output.final = data.frame(Cstorage,Mleaf,Mstem,Mroot,Sleaf)

# Plant Carbon pools for various parameter sensitivity
output.final$Date = data$Date
names(output.final) = c("Cstorage","Mleaf","Mstem","Mroot","Sleaf","Date")
melted.output = melt(output.final[,c("Mleaf","Mstem","Mroot","Cstorage","Sleaf","Date")], id.vars="Date")
melted.Cstorage = output.final[,c("Cstorage","Date")]
melted.output$Date = as.Date(melted.output$Date)
melted.Cstorage$Date = as.Date(melted.Cstorage$Date)
melted.output$Case = as.factor(q)
melted.Cstorage$Case = as.factor(q)

# Storing the summary of data, outputs, Cstorage, parameters
if (q == 1) {
  summary.output = melted.output
  summary.Cstorage = melted.Cstorage
}
if (q > 1) {
  summary.output = rbind(summary.output,melted.output)
  summary.Cstorage = rbind(summary.Cstorage,melted.Cstorage)
}

############----------------------------------------
# Take the parameters Y, k, af, as, ar, sf for 5L pot
q=5 # Case 5
param = subset(param.summary,(volume.group %in% 1))
keeps = c("Date", "variable", "Parameter")
param = param[ , keeps, drop = FALSE]
param.casted = dcast( param , Date ~ variable )

# Calculating model outputs
Mleaf = Mstem = Mroot = c()
Mleaf[1] <- data$Mleaf[1]
Mstem[1] <- data$Mstem[1]
Mroot[1] <- data$Mroot[1]

Cday=data$Cday; Rd=data$Rd;
k=param.casted$k; Y=param.casted$Y; af=param.casted$af; as=param.casted$as; ar=param.casted$ar; sf=param.casted$sf

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

for (i in 2:length(Cday)) {
  M[i-1] <- b * LA[i-1] + intercept
  GPP[i-1] <- LA[i-1] * Cday[i-1] * M[i-1] # calculate total daily C gain with self shading
  
  Cstorage[i] <- Cstorage[i-1] + GPP[i-1] - Rd[i-1]*(Mleaf[i-1] + Mroot[i-1] + Mstem[i-1]) - k*Cstorage[i-1]
  Sleaf[i] <- Cstorage[i] * 0.75 # 75% of storage goes to leaf (Duan's experiment)
  Sstem[i] <- Cstorage[i] * 0.16 # 16% of storage goes to stem (Duan's experiment)
  Sroot[i] <- Cstorage[i] * 0.09 # 9% of storage goes to root (Duan's experiment)
  
  Cleaf[i] <- Cleaf[i-1] + k*Cstorage[i-1]*af*(1-Y) - sf*Cleaf[i-1]
  Cstem[i] <- Cstem[i-1] + k*Cstorage[i-1]*as*(1-Y)
  Croot[i] <- Croot[i-1] + k*Cstorage[i-1]*(1-af-as)*(1-Y)
  
  Mleaf[i] <- Cleaf[i] + Sleaf[i]
  Mstem[i] <- Cstem[i] + Sstem[i]
  Mroot[i] <- Croot[i] + Sroot[i]
  
  # Leaf area (t) = Leaf area (T) * Leaf count (t) / Leaf count (T); t = time, T = time of harvest
  LA[i] <- leaf.data$final_LA * Cleaf[i] / leaf.data$final_LM
}
output.final = data.frame(Cstorage,Mleaf,Mstem,Mroot,Sleaf)

# Plant Carbon pools for various parameter sensitivity
output.final$Date = data$Date
names(output.final) = c("Cstorage","Mleaf","Mstem","Mroot","Sleaf","Date")
melted.output = melt(output.final[,c("Mleaf","Mstem","Mroot","Cstorage","Sleaf","Date")], id.vars="Date")
melted.Cstorage = output.final[,c("Cstorage","Date")]
melted.output$Date = as.Date(melted.output$Date)
melted.Cstorage$Date = as.Date(melted.Cstorage$Date)
melted.output$Case = as.factor(q)
melted.Cstorage$Case = as.factor(q)

# Storing the summary of data, outputs, Cstorage, parameters
if (q == 1) {
  summary.output = melted.output
  summary.Cstorage = melted.Cstorage
}
if (q > 1) {
  summary.output = rbind(summary.output,melted.output)
  summary.Cstorage = rbind(summary.Cstorage,melted.Cstorage)
}


#####################
# Plot the C pools for various test cases with parameter shifted from free to potted seedlings
summary.output.Mleaf = subset(summary.output,(variable %in% "Mleaf"))

p1 = ggplot() +
  geom_line(data = summary.output.Mleaf, aes(x = Date, y = value, group = interaction(variable,Case), colour=Case, linetype=variable)) +
  geom_point(size=2) +
  ylab("C pools (g C)") +
  ggtitle(paste("Leaf boimass")) +
  # labs(colour="C pools", linetype=as.character(q), shape=as.character(var[p])) +
  theme_bw() +
  annotate("text", x = mean(summary.output.Mleaf$Date), y = mean(summary.output.Mleaf$value), size = 3,
           label = paste("case 1 = All Parameters for free seedling",
                         "\ncase 2 = k,Y,sf for free seedling and af,as,ar for potted seedlings",
                         "\ncase 3 = k,sf for free seedling and Y,af,as,ar for potted seedlings",
                         "\ncase 4 = k for free seedling and Y,af,as,ar,sf for potted seedlings",
                         "\ncase 5 = All Parameters for potted seedlings")) +
  theme(plot.title = element_text(size = 12, face = "bold")) +
  theme(legend.title = element_text(colour="chocolate", size=12, face="bold")) +
  theme(axis.title.x = element_text(size = 12, vjust=-.2)) +
  theme(axis.title.y = element_text(size = 12, vjust=0.3))
p1
ggsave(p1,filename=paste("output/figures/Mleaf_pools.png",sep=""))

summary.output.Mstem = subset(summary.output,(variable %in% "Mstem"))

p2 = ggplot() +
  geom_line(data = summary.output.Mstem, aes(x = Date, y = value, group = interaction(variable,Case), colour=Case, linetype=variable)) +
  geom_point(size=2) +
  ylab("C pools (g C)") +
  ggtitle(paste("Stem boimass")) +
  # labs(colour="C pools", linetype=as.character(q), shape=as.character(var[p])) +
  theme_bw() +
  annotate("text", x = mean(summary.output.Mstem$Date), y = mean(summary.output.Mstem$value), size = 3,
           label = paste("case 1 = All Parameters for free seedling",
                         "\ncase 2 = k,Y,sf for free seedling and af,as,ar for potted seedlings",
                         "\ncase 3 = k,sf for free seedling and Y,af,as,ar for potted seedlings",
                         "\ncase 4 = k for free seedling and Y,af,as,ar,sf for potted seedlings",
                         "\ncase 5 = All Parameters for potted seedlings")) +
  theme(plot.title = element_text(size = 12, face = "bold")) +
  theme(legend.title = element_text(colour="chocolate", size=12, face="bold")) +
  theme(axis.title.x = element_text(size = 12, vjust=-.2)) +
  theme(axis.title.y = element_text(size = 12, vjust=0.3))
p2
ggsave(p2,filename=paste("output/figures/Mstem_pools.png",sep=""))

summary.output.Mroot = subset(summary.output,(variable %in% "Mroot"))

p3 = ggplot() +
  geom_line(data = summary.output.Mroot, aes(x = Date, y = value, group = interaction(variable,Case), colour=Case, linetype=variable)) +
  geom_point(size=2) +
  ylab("C pools (g C)") +
  ggtitle(paste("Root boimass")) +
  # labs(colour="C pools", linetype=as.character(q), shape=as.character(var[p])) +
  theme_bw() +
  annotate("text", x = mean(summary.output.Mroot$Date), y = mean(summary.output.Mroot$value), size = 3,
           label = paste("case 1 = All Parameters for free seedling",
                         "\ncase 2 = k,Y,sf for free seedling and af,as,ar for potted seedlings",
                         "\ncase 3 = k,sf for free seedling and Y,af,as,ar for potted seedlings",
                         "\ncase 4 = k for free seedling and Y,af,as,ar,sf for potted seedlings",
                         "\ncase 5 = All Parameters for potted seedlings")) +
  theme(plot.title = element_text(size = 12, face = "bold")) +
  theme(legend.title = element_text(colour="chocolate", size=12, face="bold")) +
  theme(axis.title.x = element_text(size = 12, vjust=-.2)) +
  theme(axis.title.y = element_text(size = 12, vjust=0.3))
p3
ggsave(p3,filename=paste("output/figures/Mroot_pools.png",sep=""))
