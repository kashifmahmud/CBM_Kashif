# R code to import and process Asat data, and then check the photosynthesis variation over time
rm(list=ls())

vol = c(5,10,15,20,25,35,1000)
plots = c(1,2,3,4,5,6,7)
Asat_raw <- read.csv("https://raw.githubusercontent.com/kashifmahmud/EucPVE/master/calculated%20data/Asat.csv")
keeps <- c("Date", "Photo", "volume", "plot")
Asat = Asat_raw[ , keeps, drop = FALSE]
Asat$Date = as.Date(Asat$Date)
Asat$volume = as.factor(Asat$volume)
Asat$plot = as.factor(Asat$plot)

# # plot_summary.xlsx - volume designations by plot and pot
# plot.summary = read.csv("https://raw.githubusercontent.com/kashifmahmud/EucPVE/master/raw%20data/plot_summary.csv")

# plot Asat for various treatments over time
for (i in 1:length(vol)) {
  png(file = paste("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/output/figures/Asat_vol_",vol[i],".png",sep=""))
  par(mfrow=c(3,3),oma = c(0, 0, 2, 0))
  for (j in 1:length(plots)) {
  Asat.set = subset(Asat,(volume %in% vol[i] & plot %in% plots[j]))
  plot(Asat.set$Date, Asat.set$Photo, main=paste("Asat for plot",plots[j]), xlab="Days", ylab="Asat (µmol CO2 m-2 s-1)")
  abline(lm(Photo ~ Date, data = Asat.set))
  }
  title(main = paste("Asat for volume",vol[i],"L"), outer=TRUE, cex = 1.5)
  dev.off()
}



