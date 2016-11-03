# Carbon balance model 
# Developed by Kashif Mahmud and Belinda Medlyn (October 2016)
# k.mahmud@westernsydney.edu.au

###### R script to import and process the raw HIE pot experiment data to model the carbon pools and fluxes using MCMC
# Set working directory for saving figures
setwd("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif")

# This script cleans the workspace, loads necessary Rfunctions and packages
source("load_packages_functions_CBM.R")

# Calculate Cleaf pool
# inputs:
number_of_pots = 7; # number of soil manipulation tests with volumes of 5, 10, 15, 20, 25, 35, Free(1000) Litres
number_of_plots = 7; # number of replicas

# Court's data files - Soil Manipulation test
# plot_summary.xlsx - volume designations by plot and pot
plot.summary = read.csv("rawdata/plot_summary.csv")
head(plot.summary)
# Find different soil volume pots from the data set
vols = unique(plot.summary$volume, fromLast = TRUE)
vols = vols[order(vols)]

# leaf_data.xlsx - Leaf mass and SLA data (fortnightly measurements)
leaf.data = read.csv("rawdata/leaf_data.csv")
leaf.data$Date = as.Date(leaf.data$Date, format = "%d/%m/%Y")
head(leaf.data)
class(leaf.data$Date)

# leaf_count.xlsx - Leaf count data (weekly measurements)
leaf.count = read.csv("rawdata/leaf_count_1.csv")
names(leaf.count)[1:2] <- c("plot", "pot")
lc <- merge(plot.summary,leaf.count,by=c("pot","plot"))
lc.read = read.csv("rawdata/leaf_count_1.csv", header = F)
lc.final <- as.data.frame(t(lc.read[1,3:ncol(leaf.count)]))
names(lc.final)[1] <- "Date"
lc.final$Date = as.Date(lc.final$Date, format = "%m/%d/%Y")
dimnames(lc.final)[[1]] <- 1:nrow(lc.final)
lc.SE.final = lc.final

################### Count total leaf numbers for various soil manipulation tests (Weekly data) with SE
for(i in 1:length(vols)) { #-- Create objects  'lc.idn.1', 'lc.idn.2', ... 'lc.idn.7' --
  lc.idn = subset(lc,volume==vols[i]) 
  lc.idn[nrow(lc.idn)+1, ] = colMeans(lc.idn, na.rm = TRUE) # R8 = Average of leaf counts
  lc.idn[nrow(lc.idn)+1, ] = (apply(lc.idn[1:ncol(lc.idn)], 2, sd, na.rm = TRUE))/(ncol(lc.idn)-3)^0.5 # R9 = Standard deviation of leaf counts
  # lc.idn[nrow(lc.idn)+1, ] = apply(lc.idn[1:ncol(lc.idn)], 2, sd, na.rm = TRUE) # R9 = Standard deviation of leaf counts
  dimnames(lc.idn)[[1]] <- c(1:(nrow(lc.idn)-2), "Mean", "SE")
  # lc.idn[] <- round(lc.idn)
  nam1 <- paste("volume_", vols[i], sep = "")
  nam2 <- paste("SE_volume_", vols[i], sep = "")
  lc.final[,i+1] <- as.data.frame(t(lc.idn[nrow(lc.idn)-1, 4:ncol(lc.idn)]))
  lc.SE.final[,i+1] = as.data.frame(t(lc.idn[nrow(lc.idn), 4:ncol(lc.idn)]))
  names(lc.final)[(i+1)] <- nam1
  names(lc.SE.final)[(i+1)] <- nam2
}

################### Interpolate leaf counts for daily data from weekly measurements
lc.melt <- melt(lc.final, id.vars = "Date")
names(lc.melt)[2:3] <- c("volume", "Leaf_Count")

# Perform Cubic Spline interpolation to get daily leaf count
lc.daily = data.frame(seq(lc.melt$Date[1], lc.melt$Date[nrow(lc.melt)], 1))
names(lc.daily)[1] <- "Date"
number_of_day = nrow(lc.daily)
for(i in 1:length(vols)) {
  func = splinefun(x=lc.final$Date, y=lc.final[ ,i+1], method="fmm",  ties = mean)
  lc.daily[ ,i+1] = round(func(lc.daily$Date))
}
names(lc.daily)[2:ncol(lc.daily)] <- c("volume_5", "volume_10", "volume_15", "volume_20", "volume_25", "volume_35", "volume_1000")
lc.daily.melt <- melt(lc.daily, id.vars = "Date")
names(lc.daily.melt)[2:3] <- c("volume", "Leaf_Count")
# plot interpolated daily leaf counts
# png(file = "/Users/kashifmahmud/WSU/ARC_project/CBM/Results/Leaf_count_daily.png")
ggplot() +
  geom_line(data = lc.daily.melt, aes(x = Date, y = Leaf_Count, group = volume, colour=factor(volume))) + 
  xlab("Date") +
  ylab("Daily leaf count") +
  ggtitle("Daily leaf count")
# dev.off()


################### Harvested leaf mass and leaf area
harvest_data = read.csv("rawdata/harvest aboveground mass.csv")
harvest_data <- merge(plot.summary,harvest_data,by=c("pot","plot"))
for(i in 1:length(vols)) {
  hd.idn = subset(harvest_data,volume==vols[i]) 
  hd.idn[nrow(hd.idn)+1, ] = colMeans(hd.idn, na.rm = TRUE) # R8 = Average of leaf counts
  hd.idn[nrow(hd.idn)+1, ] = (apply(hd.idn, 2, sd))/(nrow(hd.idn)-1)^0.5 # R9 = Standard deviation of leaf counts
  # hd.idn[nrow(hd.idn)+1, ] = apply(hd.idn, 2, sd) # R9 = Standard deviation of leaf counts
  hd.idn$newleaf_count <- round(hd.idn$newleaf_count)
  hd.idn$leaf_count <- round(hd.idn$leaf_count)
  hd.idn$leaf_mass <- hd.idn$Leafmass.bag - hd.idn$Leaf_bag
  hd.idn$stem_mass <- hd.idn$stemmass.bag - hd.idn$stem_bag
  keeps <- c("volume", "leaf_area", "leaf_count", "leaf_mass", "stem_mass", "Croot", "Froot")
  hd = hd.idn[ , keeps, drop = FALSE]
  dimnames(hd)[[1]] <- c(1:7, "mean", "SE")
  if (i == 1) {
    hd.final <- hd[0,]
    }
  hd.final[i, ] <- hd["mean", ]
}
dimnames(hd.final)[[1]] <- c(1:7)

################### Calculate leaf area directly from harvest data
# Leaf area (t) = Leaf area (T) * Leaf count (t) / Leaf count (T); t = time, T = time of harvest
la.daily = data.frame(lc.daily$Date)
for(i in 1:length(vols)) {
  la.daily[ , i+1] = lc.daily[ , i+1] * hd.final$leaf_area[i] / hd.final$leaf_count[i]
}
names(la.daily)[1:ncol(la.daily)] <- c("Date", "volume_5", "volume_10", "volume_15", "volume_20", "volume_25", "volume_35", "volume_1000")

la.daily.melt <- melt(la.daily, id.vars = "Date")
names(la.daily.melt)[2:3] <- c("volume", "Leaf_area")
# plot interpolated daily leaf area from harvest data
# png(file = "/Users/kashifmahmud/WSU/ARC_project/CBM/Results/Leaf_area_daily.png")
ggplot() +
  geom_line(data = la.daily.melt, aes(x = Date, y = Leaf_area, group = volume, colour=factor(volume))) +
  ylab("Leaf area (cm^2)") +
  xlab("Days") +
  ggtitle("Leaf area from harvest data")
# dev.off()


################### Calculate leaf mass (leaf carbon pool) directly from harvest data
# Leaf mass (t) = Leaf mass (T) * Leaf count (t) / Leaf count (T); t = time, T = time of harvest
lm.daily = data.frame(lc.daily$Date)
for(i in 1:length(vols)) {
  lm.daily[ , i+1] = lc.daily[ , i+1] * hd.final$leaf_mass[i] / hd.final$leaf_count[i]
  # lm.daily[ , i+1] = lm.daily[ , i+1] * 0.5 # Unit connversion = gm of DM from gm of C (multiplying by 0.5)
}
names(lm.daily)[1:ncol(lm.daily)] <- c("Date", "volume_5", "volume_10", "volume_15", "volume_20", "volume_25", "volume_35", "volume_1000")

lm.daily.melt <- melt(lm.daily, id.vars = "Date")
names(lm.daily.melt)[2:3] <- c("volume", "leafmass")
# plot interpolated daily leaf mass from harvest data
# png(file = "/Users/kashifmahmud/WSU/ARC_project/CBM/Results/Leaf_mass_daily.png")
ggplot() +
  geom_line(data = lm.daily.melt, aes(x = Date, y = leafmass, group = volume, colour=factor(volume))) +
  ylab("Leaf mass (gmC)") +
  xlab("Days") +
  ggtitle("Leaf mass from harvest data")
# dev.off()

# Save the weekly Cleaf data for MCMC CBM
lm.weekly = subset(lm.daily,Date %in% lc.final$Date)
lm.weekly.melt <- melt(lm.weekly, id.vars = "Date")
names(lm.weekly.melt)[2:3] <- c("volume", "leafmass")
lm.weekly.melt$volume = as.numeric(lm.weekly.melt$volume)
for(i in length(vols):1) {
  ind = which(lm.weekly.melt$volume %in% i)
  lm.weekly.melt$volume[ind] = vols[i]
}
# lm.SE.weekly = subset(lm.SE.daily,Date %in% lc.SE.final$Date)
lc.weekly.melt <- melt(lc.final, id.vars = "Date")
names(lc.weekly.melt)[2:3] <- c("volume", "leafcount")
lc.weekly.melt$volume = as.numeric(lc.weekly.melt$volume)
for(i in length(vols):1) {
  ind = which(lc.weekly.melt$volume %in% i)
  lc.weekly.melt$volume[ind] = vols[i]
}
lc.SE.weekly.melt <- melt(lc.SE.final, id.vars = "Date")
names(lc.SE.weekly.melt)[2:3] <- c("volume", "leafcount_SE")
lc.SE.weekly.melt$volume = as.numeric(lc.SE.weekly.melt$volume)
for(i in length(vols):1) {
  ind = which(lc.SE.weekly.melt$volume %in% i)
  lc.SE.weekly.melt$volume[ind] = vols[i]
}
# leaf.data.gas = merge(sla.final,lc.daily.m,by=c("Date","volume"))
lm.weekly.melt = merge(lm.weekly.melt,lc.weekly.melt,by=c("Date","volume"))
lm.weekly.melt = merge(lm.weekly.melt,lc.SE.weekly.melt,by=c("Date","volume"))
lm.weekly.melt$leafmass_SE = lm.weekly.melt$leafmass * lm.weekly.melt$leafcount_SE / lm.weekly.melt$leafcount
write.csv(lm.weekly.melt[,c("Date","volume","leafmass","leafmass_SE")], file = "rawdata/Cleaf_weekly_data.csv", row.names = FALSE)


################### Analyse stem height diameter to estimate Stem carbon pool
# Import weekly height diameter data for 3 months (Diameter is in mm; Height is in cm)
height.dia <- read.csv("rawdata/height_diameter.csv")
library(lubridate)
height.dia$Date <- parse_date_time(height.dia$Date,"d m y")
height.dia$Date = as.Date(height.dia$Date, format = "%d/%m/%Y")
height.dia <- merge(plot.summary,height.dia,by=c("pot","plot"))
height.dia = height.dia[ ,c(-1,-2)]
height.dia = height.dia[with(height.dia, order(volume,Date)), ]
height.dia <- height.dia[c("Date", "volume", "height", "diameter")]

for(i in 1:length(vols)) {
  height.dia.idn = subset(height.dia,volume==vols[i]) 
  for(j in 1:length(unique(height.dia.idn$Date))) {
    height.dia.idn.date = subset(height.dia.idn, Date == unique(height.dia.idn$Date)[j])
    height.dia.idn.date[nrow(height.dia.idn.date)+1, 2:ncol(height.dia.idn.date)] = colMeans(height.dia.idn.date[2:ncol(height.dia.idn.date)], na.rm = TRUE) # R7 = Average of leaf data
    height.dia.idn.date[nrow(height.dia.idn.date)+1, 2:ncol(height.dia.idn.date)] = (apply(height.dia.idn.date[2:ncol(height.dia.idn.date)], 2, sd, na.rm = TRUE))/(nrow(height.dia.idn.date)-1)^0.5 # R8 = Standard deviation of leaf counts
    # height.dia.idn.date[nrow(height.dia.idn.date)+1, 2:ncol(height.dia.idn.date)] = apply(height.dia.idn.date[2:ncol(height.dia.idn.date)], 2, sd, na.rm = TRUE) # R8 = Standard deviation of leaf counts
    height.dia.idn.date$Date = height.dia.idn.date$Date[1]
    dimnames(height.dia.idn.date)[[1]] <- c(1:(nrow(height.dia.idn.date)-2), "Mean", "SE")
    if (i == 1 && j == 1) {
      height.dia.final <- height.dia.idn.date[0,]
    }
    height.dia.final[j+(i-1)*length(unique(height.dia.idn$Date)), ] <- height.dia.idn.date["Mean", ]
    height.dia.final$height_SE[j+(i-1)*length(unique(height.dia.idn$Date))] <- height.dia.idn.date["SE", 3]
    height.dia.final$dia_SE[j+(i-1)*length(unique(height.dia.idn$Date))] <- height.dia.idn.date["SE", 4]
  }
}

################### Linear regression model fitting [log(wood_mass) = b(1) + b(2)*log(dia) + b(3)*log(height)]
# Fit the model with initial data (10) and harvested data for free seedlings (7)
# Import initial seedling data
initial.data <- read.csv("rawdata/seedling_initial.csv")
initial.data[nrow(initial.data)+1, 2:ncol(initial.data)] = colMeans(initial.data[2:ncol(initial.data)], na.rm = TRUE) # R7 = Average of leaf data
initial.data[nrow(initial.data)+1, 2:ncol(initial.data)] = (apply(initial.data[2:ncol(initial.data)], 2, sd, na.rm = TRUE))/(nrow(initial.data)-1)^0.5 # R8 = Standard deviation of leaf counts
# initial.data[nrow(initial.data)+1, 2:ncol(initial.data)] = apply(initial.data[2:ncol(initial.data)], 2, sd, na.rm = TRUE) # R8 = Standard deviation of leaf counts
dimnames(initial.data)[[1]] <- c(1:(nrow(initial.data)-2), "Mean", "SE")

# Import harvested seedling data for all different treatments
end.data <- read.csv("rawdata/seedling_mass.csv")
end.data = end.data[ ,c(-1,-2)]
end.data = end.data[with(end.data, order(volume)), ]
for(i in 1:length(vols)) {
  end.data.idn = subset(end.data,volume==vols[i]) 
  end.data.idn[nrow(end.data.idn)+1, 2:ncol(end.data.idn)] = colMeans(end.data.idn[2:ncol(end.data.idn)], na.rm = TRUE) # R7 = Average of leaf data
  end.data.idn[nrow(end.data.idn)+1, 2:ncol(end.data.idn)] = (apply(end.data.idn[2:ncol(end.data.idn)], 2, sd, na.rm = TRUE))/(nrow(end.data.idn)-1)^0.5 # R8 = Standard deviation of leaf counts
  # end.data.idn[nrow(end.data.idn)+1, 2:ncol(end.data.idn)] = apply(end.data.idn[2:ncol(end.data.idn)], 2, sd, na.rm = TRUE) # R8 = Standard deviation of leaf counts
  end.data.idn$volume = end.data.idn$volume[1]
  dimnames(end.data.idn)[[1]] <- c(1:(nrow(end.data.idn)-2), "Mean", "SE")
  if (i == 1) {
    end.data.final <- end.data.idn[0,]
  }
  end.data.final[i, ] <- end.data.idn["Mean", ]
  end.data.final$coarseroot_SE[i] <- end.data.idn["SE", 2]
  end.data.final$fineroot_SE[i] <- end.data.idn["SE", 3]
  end.data.final$leafmass_SE[i] <- end.data.idn["SE", 4]
  end.data.final$stemmass_SE[i] <- end.data.idn["SE", 5]
  end.data.final$totalmass_SE[i] <- end.data.idn["SE", 6]
}

# Linear regression model fitting: log(stem_mass) = b(1) + b(2)*log(dia) + b(3)*log(height)
stemmass = c(initial.data$wood_mass[1:10], hd$stem_mass[1:7])
p = subset(height.dia.idn,Date=="2013-05-21")
height = c(initial.data$height[1:10], p$height)
diameter = c(initial.data$diameter_15[1:10], p$diameter)

model.fit = data.frame(stemmass, height, diameter)
fit <- lm(log(stemmass) ~ log(diameter) + log(height), data=model.fit)
summary(fit) # show results
coefficients(fit) # model coefficients
cat("Linear regression model fitting: log(stem_mass) = ", coefficients(fit)[1], "+", coefficients(fit)[2], 
    "* log(diameter) +", coefficients(fit)[3], "* log(height)")

# Estimate the stemmass from the fitted linear regression equation
eq = function(x,y){exp(coefficients(fit)[1] + coefficients(fit)[2] * log(x)  + coefficients(fit)[3] * log(y))}
x = height.dia.final$diameter
y = height.dia.final$height
z = height.dia.final$stemmass = eq(x,y)


# Calculate all seedling stem mass from height and diameter using the linear model and then get the SEs from the 7 replicas
height.dia$stemmass = eq(height.dia$diameter,height.dia$height)
for(i in 1:length(vols)) {
  stemmass.idn = subset(height.dia,volume==vols[i]) 
  for(j in 1:length(unique(stemmass.idn$Date))) {
    stemmass.idn.date = subset(stemmass.idn, Date == unique(stemmass.idn$Date)[j])
    stemmass.idn.date[nrow(stemmass.idn.date)+1, 2:ncol(stemmass.idn.date)] = colMeans(stemmass.idn.date[2:ncol(stemmass.idn.date)], na.rm = TRUE) # R7 = Average of leaf data
    stemmass.idn.date[nrow(stemmass.idn.date)+1, 2:ncol(stemmass.idn.date)] = (apply(stemmass.idn.date[2:ncol(stemmass.idn.date)], 2, sd, na.rm = TRUE))/(nrow(stemmass.idn.date)-1)^0.5 # R8 = Standard deviation of leaf counts
    # stemmass.idn.date[nrow(stemmass.idn.date)+1, 2:ncol(stemmass.idn.date)] = apply(stemmass.idn.date[2:ncol(stemmass.idn.date)], 2, sd, na.rm = TRUE) # R8 = Standard deviation of leaf counts
    stemmass.idn.date$Date = stemmass.idn.date$Date[1]
    dimnames(stemmass.idn.date)[[1]] <- c(1:(nrow(stemmass.idn.date)-2), "Mean", "SE")
    if (i == 1 && j == 1) {
      stemmass.final <- stemmass.idn.date[0,]
    }
    stemmass.final[j+(i-1)*length(unique(stemmass.idn$Date)), ] <- stemmass.idn.date["Mean", ]
    stemmass.final$height_SE[j+(i-1)*length(unique(stemmass.idn$Date))] <- stemmass.idn.date["SE", 3]
    stemmass.final$dia_SE[j+(i-1)*length(unique(stemmass.idn$Date))] <- stemmass.idn.date["SE", 4]
    stemmass.final$stemmass_SE[j+(i-1)*length(unique(stemmass.idn$Date))] <- stemmass.idn.date["SE", 5]
  }
}

# Save the weekly Cstem data for MCMC CBM
write.csv(stemmass.final[ , c("Date","volume","stemmass","stemmass_SE")], file = "rawdata/Cstem_weekly_data.csv", row.names = FALSE)


# Calculate SE for initial and harvest rootmass data
Croot = data.frame(Date = as.Date(c("2013-01-21","2013-05-21")), rootmass = numeric(2), rootmass_SE = numeric(2))
Croot = merge(vols,Croot)
names(Croot)[1] = "volume"
Croot = Croot[with(Croot, order(Date)), ]
end.data.final$rootmass = end.data.final$coarseroot + end.data.final$fineroot
end.data.final$rootmass_SE = end.data.final$coarseroot_SE + end.data.final$fineroot_SE

Croot$rootmass[Croot$Date=="2013-01-21"] = initial.data$root_mass[11]
Croot$rootmass_SE[Croot$Date=="2013-01-21"] = initial.data$root_mass[12]
Croot$rootmass[Croot$Date=="2013-05-21"] = end.data.final$rootmass
Croot$rootmass_SE[Croot$Date=="2013-05-21"] = end.data.final$rootmass_SE
# Save the initial and harvest rootmass data for MCMC CBM
write.csv(Croot, file = "rawdata/Croot_twice_data.csv", row.names = FALSE)


# Plotting observation and modelled data
# png(file = "/Users/kashifmahmud/WSU/ARC_project/CBM/Results/Height_Stem mass.png")
plot(model.fit$height,model.fit$stemmass,col="red",main="Height vs Stem mass", pch=15)
lines(y,z,type="p",xlab="Height (cm)", ylab="Stem mass (g)",col="green", pch=20)
legend('topleft', c("Measurements", "Modelled data"), lty=1, col=c('red','green'), bty='n', cex=0.75)
# dev.off()
# png(file = "/Users/kashifmahmud/WSU/ARC_project/CBM/Results/Diameter_Stem mass.png")
plot(model.fit$diameter,model.fit$stemmass, col="red",main="Diameter vs Stem mass")
lines(x,z,type="p",col="green", xlab="Diameter (cm)", ylab="Stem mass (g)")
legend('topleft', c("Measurements", "Modelled data"), lty=1, col=c('red','green'), bty='n', cex=0.75)
# dev.off()


################### Calculate daily GPP
##read data
finalmass = read.csv("rawdata/harvest_mass_means.csv")

la.daily.m = la.daily.melt # Leaf area calculated from harvested data
la.daily.m$volume = as.numeric(la.daily.m$volume)
for(i in length(vols):1) {
  ind = which(la.daily.m$volume %in% i)
  la.daily.m$volume[ind] = vols[i]
}
la.daily.m$Leaf_area = la.daily.m$Leaf_area / (100*100) # unit conversion from cm2 to m2

Cday_gross <- read.csv("rawdata/cday_120_clean_gross.csv")
Cday_gross$Date <- as.Date(Cday_gross$Date)

##Cday with leaf area needs self shading (use slope intercept, from 5-free vol)
sigma <- read.csv("rawdata/M_leafarea_model.csv")

##function to generate total plant daily C gain-----------------------------------------------------------------------------
modelledC_func <- function(leafarea, shading, Cday){ #leafarea dfr, self shading dfr, and modelled C gain (gm2)
  dailyCnet <- merge(leafarea,Cday)
  # this needs to include self shadeing (M as a linear function of leaf area)
  dailyCnet <- merge(dailyCnet, shading[, c(2,3,5)], by="volume")
  dailyCnet$M <- with(dailyCnet, b*Leaf_area+intercept)
  #calculate total daily C gain with self shading
  dailyCnet$tdcg <- with(dailyCnet, Leaf_area * carbon_day * M)
  print("successfully calculated total daily carbon with modelled Cgain, leaf area and self shading")
  
  dailyCnet$volume <- as.factor(dailyCnet$volume)
  return(dailyCnet)
}

##calculate daily gross C gain
dailyCgross <- modelledC_func(la.daily.m, sigma, Cday_gross)
names(dailyCgross)[8] <- "tdc_gross"
names(dailyCgross)[4] <- "cday_gross"

# dailyC <- merge(dailyCgross[,c(1:2,4,8)], dailyCnet[,c(1:2,4,8)], by=c("Date", "volume")) # Unit = gC d-1
write.csv(dailyCgross[ , c("Date","volume","tdc_gross")], "rawdata/GPP.csv", row.names = FALSE)

ggplot(data = dailyCgross, aes(x = Date, y = tdc_gross, group = volume, colour=factor(volume))) + 
  geom_point(size=1) +
  xlab("Date") +
  ylab("Gross total daily C (gC)") +
  ggtitle("Gross total daily C")


#################### Rdark prediction through time using rdarkq10 equation by volume
#import site weather data, rename variables, format date stuff
eucpve_met <- read.csv("rawdata/eucpve_met.csv")
eucpve_met <- eucpve_met[ , c(1,4)]
names(eucpve_met)[2] <- "temp"
eucpve_met$Date <- as.Date(eucpve_met$DateTime15)

#need to turn the datetime 15 into hms
#install.packages("lubridate")
eucpve_met$DateTime15 <- ymd_hms(eucpve_met$DateTime15)
eucpve_met$time15 <- format(eucpve_met$DateTime15, format='%H:%M:%S')

#subset by Date range of experiment
eucpve_met1 <- subset(eucpve_met[, 2:4], Date  >= "2013-01-21" & Date  <= "2013-05-21")

#Rdark Q10 equations by volume
rdarkq10 <- read.csv("rawdata/rdarkq10.csv")

### A_model, use temperature data from site to generate daily total leaf dark respiration to enter into the exponential Rd model
volume <- data.frame(vols)
A_model <- merge(eucpve_met1, volume)
names(A_model)[4] <- "volume"
#need to calculate Rdark through time using rdarkq10 equation by volume
A_model <- merge(A_model, rdarkq10[,1:2], by="volume")
q25_drake <- 1.86

# refit RD (gC m-2 leaf d-1)
hd.final$sla_harvest = hd.final$leaf_area / hd.final$leaf_mass / 100^2 # calculate harvested SLA
A_model <- merge(A_model, hd.final[,c(1,8)], by="volume")
A_model$Rd_daily <- with(A_model, rd12.3 * q25_drake^((temp-12.3)/10) * sla_harvest) # unit (micromol CO2 g-1 plant d-1)
with(A_model, plot(temp,Rd_daily, col=volume))


Rd <- summaryBy(Rd_daily ~ Date+volume, data=A_model, FUN=sum, keep.names=TRUE ) # Sum of all same day Rd
Rd$Rd_daily = Rd$Rd_daily / (24*4) # Average of all 96 data from a day
Rd$Rd_daily = Rd$Rd_daily * (12/44)  # unit conversion from micromol CO2 to gC
write.csv(Rd, "rawdata/Rd.csv", row.names=FALSE)


################### Find plant storage (tnc, fortnightly data) for corresponding Dates (from Court's leaf_data file: represents Gas measurement campaign)
keeps <- c("Date", "volume", "starch_mgperg", "sugars_mgperg")
tnc.gas = leaf.data[ , keeps, drop = FALSE]
names(tnc.gas)[1:2] <- c("Date", "volume")
tnc.gas$tnc = tnc.gas$starch_mgperg + tnc.gas$sugars_mgperg
keeps <- c("Date", "volume", "tnc") # Unit = mg g-1leaf
tnc.gas = tnc.gas[ , keeps, drop = FALSE]

for(i in 1:length(vols)) {
  tnc.idn = subset(tnc.gas,volume==vols[i]) 
  for(j in 1:length(unique(tnc.idn$Date))) {
    tnc.idn.date = subset(tnc.idn, Date == unique(tnc.idn$Date)[j])
    tnc.idn.date[nrow(tnc.idn.date)+1, 2:ncol(tnc.idn.date)] = colMeans(tnc.idn.date[2:ncol(tnc.idn.date)], na.rm = TRUE) # R7 = Average of tnc
    tnc.idn.date[nrow(tnc.idn.date)+1, 2:ncol(tnc.idn.date)] = (apply(tnc.idn.date[2:ncol(tnc.idn.date)], 2, sd))/(nrow(tnc.idn.date)-1)^0.5 # R8 = Standard deviation of tnc
    # tnc.idn.date[nrow(tnc.idn.date)+1, 2:ncol(tnc.idn.date)] = apply(tnc.idn.date[2:ncol(tnc.idn.date)], 2, sd) # R8 = Standard deviation of tnc
    tnc.idn.date$Date = tnc.idn.date[1,1]
    dimnames(tnc.idn.date)[[1]] <- c(1:(nrow(tnc.idn.date)-2), "Mean", "SE")
    if (i == 1 && j == 1) {
      tnc.final <- tnc.idn.date[0,]
    }
    tnc.final[j+(i-1)*length(unique(tnc.idn$Date)), ] <- tnc.idn.date["Mean", ]
    tnc.final$tnc_SE[j+(i-1)*length(unique(tnc.idn$Date))] <- tnc.idn.date["SE", 3]
  }
}
# Unit conversion from (mg g-1leaf) to  (g plant-1)
lm.daily.m = lm.daily.melt # Leaf area calculated from harvested data
lm.daily.m$volume = as.numeric(lm.daily.m$volume)
for(i in length(vols):1) {
  ind = which(lm.daily.m$volume %in% i)
  lm.daily.m$volume[ind] = vols[i]
}
lm.daily.m$Date = as.Date(lm.daily.m$Date)
lm.daily.gas = lm.daily.m[lm.daily.m$Date %in% as.Date(c(unique(tnc.final$Date))), ]
tnc.final$tnc = tnc.final$tnc * lm.daily.gas$leafmass / 1000 # Unit = g plant-1
tnc.final$tnc_SE = tnc.final$tnc_SE * lm.daily.gas$leafmass / 1000 # Unit = g plant-1
write.csv(tnc.final, file = "rawdata/tnc_fortnightly_data.csv", row.names = FALSE)



