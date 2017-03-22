# R code to import and process Asat data, and then plot the GPP, Asat, N content, TNC, Starch content variation over time
# Plot weekly cumulative GPP and Leaf mass for various treatments

rm(list=ls())

vols = c(5,10,15,20,25,35,1000)
# plots = c(1,2,3,4,5,6,7)
# Import and process the Asat data
Asat_raw <- read.csv("https://raw.githubusercontent.com/kashifmahmud/EucPVE/master/calculated%20data/Asat.csv")
keeps <- c("Date", "volume", "Photo") # Photo unit = µmol CO2 m-2 s-1
Asat = Asat_raw[ , keeps, drop = FALSE]
Asat$Date = as.Date(Asat$Date)
Asat$volume = as.numeric(Asat$volume)
names(Asat) = c("Date", "volume", "Asat")

for(i in 1:length(vols)) {
  Asat.idn = subset(Asat,volume==vols[i]) 
  for(j in 1:length(unique(Asat.idn$Date))) {
    Asat.idn.date = subset(Asat.idn, Date == unique(Asat.idn$Date)[j])
    Asat.idn.date[nrow(Asat.idn.date)+1, 2:ncol(Asat.idn.date)] = colMeans(Asat.idn.date[2:ncol(Asat.idn.date)], na.rm = TRUE) # R7 = Average of Asat
    Asat.idn.date[nrow(Asat.idn.date)+1, 2:ncol(Asat.idn.date)] = (apply(Asat.idn.date[2:ncol(Asat.idn.date)], 2, sd))/(nrow(Asat.idn.date)-1)^0.5 # R8 = Standard deviation of Asat
    Asat.idn.date$Date = Asat.idn.date[1,1]
    dimnames(Asat.idn.date)[[1]] <- c(1:(nrow(Asat.idn.date)-2), "Mean", "SE")
    if (i == 1 && j == 1) {
      Asat.final <- Asat.idn.date[0,]
    }
    Asat.final[j+(i-1)*length(unique(Asat.idn$Date)), ] <- Asat.idn.date["Mean", ]
    Asat.final$Asat_SE[j+(i-1)*length(unique(Asat.idn$Date))] <- Asat.idn.date["SE", 3]
  }
}

# Import and process the N content data
nc_raw <- read.csv("https://raw.githubusercontent.com/kashifmahmud/EucPVE/master/calculated%20data/leaf%20N%20content.csv")
keeps <- c("Date", "volume", "Narea") # Narea unit = µmol N m-2
nc = nc_raw[ , keeps, drop = FALSE]
nc$Date = as.Date(nc$Date)
nc$volume = as.numeric(nc$volume)
names(nc) = c("Date", "volume", "nc")

for(i in 1:length(vols)) {
  nc.idn = subset(nc,volume==vols[i]) 
  for(j in 1:length(unique(nc.idn$Date))) {
    nc.idn.date = subset(nc.idn, Date == unique(nc.idn$Date)[j])
    nc.idn.date[nrow(nc.idn.date)+1, 2:ncol(nc.idn.date)] = colMeans(nc.idn.date[2:ncol(nc.idn.date)], na.rm = TRUE) # R7 = Average of nc
    nc.idn.date[nrow(nc.idn.date)+1, 2:ncol(nc.idn.date)] = (apply(nc.idn.date[2:ncol(nc.idn.date)], 2, sd))/(nrow(nc.idn.date)-1)^0.5 # R8 = Standard deviation of nc
    nc.idn.date$Date = nc.idn.date[1,1]
    dimnames(nc.idn.date)[[1]] <- c(1:(nrow(nc.idn.date)-2), "Mean", "SE")
    if (i == 1 && j == 1) {
      nc.final <- nc.idn.date[0,]
    }
    nc.final[j+(i-1)*length(unique(nc.idn$Date)), ] <- nc.idn.date["Mean", ]
    nc.final$nc_SE[j+(i-1)*length(unique(nc.idn$Date))] <- nc.idn.date["SE", 3]
  }
}

# Import final TNC data from initial data processing
tnc.final <- read.csv("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/rawdata/tnc_fortnightly_data.csv")
# tnc unit = g plant-1
tnc.final$Date = as.Date(tnc.final$Date)
tnc.final$volume = as.numeric(tnc.final$volume)

# data = merge(Asat.final,tnc.final,by=c("Date","volume"))
# data = merge(data,nc.final,by=c("Date","volume"))

################### Find plant starch content (fortnightly data) for corresponding Dates (from Court's leaf_data file: represents Gas measurement campaign)
leaf.data <- read.csv("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/rawdata/leaf_data.csv")
leaf.data$Date = as.Date(leaf.data$Date, format = "%d/%m/%Y")
keeps <- c("Date", "volume", "starch_mgperg")
starch.gas = leaf.data[ , keeps, drop = FALSE]
names(starch.gas)[1:3] <- c("Date", "volume", "starch")
starch.gas$Date = as.Date(starch.gas$Date)
starch.gas$volume = as.numeric(starch.gas$volume)

for(i in 1:length(vols)) {
  starch.idn = subset(starch.gas,volume==vols[i]) 
  for(j in 1:length(unique(starch.idn$Date))) {
    starch.idn.date = subset(starch.idn, Date == unique(starch.idn$Date)[j])
    starch.idn.date[nrow(starch.idn.date)+1, 2:ncol(starch.idn.date)] = colMeans(starch.idn.date[2:ncol(starch.idn.date)], na.rm = TRUE) # R7 = Average of starch
    starch.idn.date[nrow(starch.idn.date)+1, 2:ncol(starch.idn.date)] = (apply(starch.idn.date[2:ncol(starch.idn.date)], 2, sd))/(nrow(starch.idn.date)-1)^0.5 # R8 = Standard error of starch
    starch.idn.date$Date = starch.idn.date[1,1]
    dimnames(starch.idn.date)[[1]] <- c(1:(nrow(starch.idn.date)-2), "Mean", "SE")
    if (i == 1 && j == 1) {
      starch.final <- starch.idn.date[0,]
    }
    starch.final[j+(i-1)*length(unique(starch.idn$Date)), ] <- starch.idn.date["Mean", ]
    starch.final$starch_SE[j+(i-1)*length(unique(starch.idn$Date))] <- starch.idn.date["SE", 3]
  }
}
# Unit conversion from (mg g-1leaf) to  (g plant-1)
leafmass.data <- read.csv("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/rawdata/leafmass_tnc_data.csv")
leafmass.data$volume = as.numeric(leafmass.data$volume)
leafmass.data$Date = as.Date(leafmass.data$Date)
starch.final = merge(starch.final,leafmass.data, by=c("Date","volume"))
starch.final$starch = starch.final$starch * starch.final$leafmass / 1000 # Unit = g plant-1
starch.final$starch_SE = starch.final$starch_SE * starch.final$leafmass / 1000 # Unit = g plant-1
# write.csv(starch.final, file = "rawdata/starch_fortnightly_data.csv", row.names = FALSE)

# Import and process the GPP data
GPP.data <- read.csv("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/rawdata/GPP.csv")
GPP.data$Date = as.Date(GPP.data$Date)
GPP.final = GPP.data[GPP.data$Date %in% as.Date(c(unique(starch.final$Date))), ]


# plot Asat for various treatments over time
png(file = paste("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/output/figures/Asat_over_time.png"))
par(mfrow=c(3,3),oma = c(0, 0, 2, 0))
for (i in 1:length(vols)) {
  Asat.set = subset(Asat.final,(volume %in% vols[i]))
  plot(Asat.set$Date, Asat.set$Asat, main=paste("Asat for volume",vols[i]), 
       ylim=c(min(Asat.final$Asat), max(Asat.final$Asat)), xlab="Days", ylab="Asat (µmol CO2 m-2 s-1)")
  abline(lm(Asat ~ Date, data = Asat.set))
}
title(main = paste("Asat for different treatments"), outer=TRUE, cex = 1.5)
dev.off()


# plot GPP for various treatments over time
png(file = paste("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/output/figures/GPP_over_time.png"))
par(mfrow=c(3,3),oma = c(0, 0, 2, 0))
for (i in 1:length(vols)) {
  GPP.set = subset(GPP.final,(volume %in% vols[i]))
  plot(GPP.set$Date, GPP.set$tdc_gross, main=paste("GPP for volume",vols[i]), 
       ylim=c(min(GPP.final$tdc_gross), max(GPP.final$tdc_gross)), xlab="Days", ylab="GPP (g C plant-1)")
  abline(lm(tdc_gross ~ Date, data = GPP.set))
}
title(main = paste("GPP for different treatments"), outer=TRUE, cex = 1.5)
dev.off()


# plot TNC for various treatments over time
png(file = paste("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/output/figures/TNC_over_time.png"))
par(mfrow=c(3,3),oma = c(0, 0, 2, 0))
for (i in 1:length(vols)) {
  tnc.set = subset(tnc.final,(volume %in% vols[i]))
  plot(tnc.set$Date, tnc.set$tnc, main=paste("TNC for volume",vols[i]), 
       ylim=c(min(tnc.final$tnc), max(tnc.final$tnc)), xlab="Days", ylab="TNC (g plant-1)")
  abline(lm(tnc ~ Date, data = tnc.set))
}
title(main = paste("TNC for different treatments"), outer=TRUE, cex = 1.5)
dev.off()


# plot Starch for various treatments over time
png(file = paste("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/output/figures/Starch_over_time.png"))
par(mfrow=c(3,3),oma = c(0, 0, 2, 0))
for (i in 1:length(vols)) {
  starch.set = subset(starch.final,(volume %in% vols[i]))
  plot(starch.set$Date, starch.set$starch, main=paste("Starch for volume",vols[i]), 
       ylim=c(min(starch.final$starch), max(starch.final$starch)), xlab="Days", ylab="Starch (g plant-1)")
  abline(lm(starch ~ Date, data = starch.set))
}
title(main = paste("Starch for different treatments"), outer=TRUE, cex = 1.5)
dev.off()


# plot N content for various treatments over time
png(file = paste("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/output/figures/N_content_over_time.png"))
par(mfrow=c(3,3),oma = c(0, 0, 2, 0))
for (i in 1:length(vols)) {
  nc.set = subset(nc.final,(volume %in% vols[i]))
  plot(nc.set$Date, nc.set$nc, main=paste("N content for volume",vols[i]), 
       ylim=c(min(nc.final$nc), max(nc.final$nc)), xlab="Days", ylab="N content (µmol m-2)")
  abline(lm(nc ~ Date, data = nc.set))
}
title(main = paste("N content for different treatments"), outer=TRUE, cex = 1.5)
dev.off()

###################################
# Import daily GPP, Rd and weekly Cleaf data
GPP.data.raw = read.csv("rawdata/GPP.csv") # Units gC d-1
# Rd.data.raw = read.csv("rawdata/Rd.csv") # Units g C g-1 plant d-1
names(GPP.data.raw) = c("Date", "volume", "gpp")
# GPP.data.raw = merge(GPP.data.raw,Rd.data.raw,by=c("Date","volume"))
GPP.data.raw$Date = as.Date(GPP.data.raw$Date)
Mleaf.data.raw = read.csv("rawdata/Cleaf_weekly_data.csv") # Units gC d-1
Mleaf.data.raw = Mleaf.data.raw[with(Mleaf.data.raw, order(Date)), ]
Mleaf.data.raw$Date = as.Date(Mleaf.data.raw$Date)

# Calculate weekly accumulated leaf mass for all treatments
for(i in 1:length(vols)) {
  Mleaf.data.idn = subset(Mleaf.data.raw,volume==vols[i]) 
  Mleaf.data.set = Mleaf.data.idn[2:nrow(Mleaf.data.idn),c("Date","volume")]
  Mleaf.data.set$weekly.LM = diff(Mleaf.data.idn$leafmass)
  if (i == 1) {
    Mleaf.data.final <- Mleaf.data.set
  }
  if (i > 1) {
    Mleaf.data.final <- rbind(Mleaf.data.final,Mleaf.data.set)
  }
}

Mleaf.data.idn$Date = as.Date(Mleaf.data.idn$Date)
GPP.data.raw$Date = as.Date(GPP.data.raw$Date)
# Calculate weekly accumulated GPP for all treatments
GPP.idn.set = GPP.idn[1,]
for(i in 1:length(vols)) {
  GPP.idn = subset(GPP.data.raw,volume==vols[i]) 
  GPP.idn = GPP.idn[with(GPP.idn, order(Date)), ]
  for(j in 1:length(unique(nc.idn$Date))) {
    GPP.idn$gpp.cum = cumsum(GPP.idn$gpp)
    for(k in 1:nrow(Mleaf.data.idn)) {
      GPP.idn.set[k,] = with(GPP.idn, GPP.idn[(Date == Mleaf.data.idn$Date[k]), ])
    }
    GPP.set = GPP.idn.set[2:nrow(GPP.idn.set),c("Date","volume")]
    GPP.set$weekly.gpp = diff(GPP.idn.set$gpp.cum)
  }
  if (i == 1) {
    GPP.final <- GPP.set
  }
  if (i > 1) {
    GPP.final <- rbind(GPP.final,GPP.set)
  }
}

data = merge(GPP.final,Mleaf.data.final,by=c("Date","volume"))

# plot interpolated daily leaf mass from harvest data
for (i in 1:length(vols)) {
  data.idn = subset(data,(volume %in% vols[i]))
  data.melt = melt(data.idn[c("Date","weekly.gpp","weekly.LM")], id.vars = "Date")
  p1 = ggplot() +
    geom_point(data = data.melt, aes(x = Date, y = value, group = variable, colour=factor(variable))) +
    geom_line(data = data.melt, aes(x = Date, y = value, group = variable, colour=factor(variable))) +
    ylab("GPP and Leaf mass (g C)") +
    xlab("Date") +
    ggtitle(paste("Volume",vols[i],"L"))
  p1
  ggsave(p1,filename=paste("output/figures/GPP and Leaf mass_vol_",vols[i],".png",sep=""))
}

# # plot weekly cumulative GPP and Leaf mass for various treatments
# png(file = paste("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/output/figures/GPP_vs_LM_weekly.png"))
# par(mfrow=c(1,1))
#     for (i in 1:length(vols)) {
#   data.idn = subset(data,(volume %in% vols[i]))
#   plot(data.idn$Date,data.idn$weekly.gpp,col="red",main=paste("Volume",vols[i]))
#   lines(data.idn$Date,data.idn$weekly.LM,type="p",xlab="Date", ylab="GPP and Leaf Mass (g C)",col="green")
#   legend('topleft', c("Measurements", "Modelled data"), lty=1, col=c('red','green'), bty='n', cex=0.75)
#   
#   
#   plot(data.idn$Date, data.idn$weekly.gpp, col="red", main=paste("Volume",vols[i]), 
#        xlab="Days", ylab="GPP and Leaf Mass")
#   lines(data.idn$Date, data.idn$weekly.LM, ylim=c(min(data.idn$weekly.LM), max(data.idn$weekly.gpp)))
# }
# title(main = paste("Weekly cumulative GPP and Leaf mass for different treatments"), outer=TRUE, cex = 1.5)
# dev.off()
