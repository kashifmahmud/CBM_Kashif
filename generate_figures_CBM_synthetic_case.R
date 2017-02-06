# Carbon balance model 
# Developed by Kashif Mahmud and Belinda Medlyn (November 2016)
# k.mahmud@westernsydney.edu.au

# This script creates the figures and saves those
##############################

# Set working directory for saving figures
setwd("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/output/figures/summary/Cpools")

# Plot modelled vs measured data ("Mleaf","Mstem","Mroot","Sleaf") against "volume" and "Total No of param"
meas = as.factor(c("Mleaf","Mstem","Mroot","Sleaf"))
res = as.factor(c("Mleaf.modelled","Mstem.modelled","Mroot.modelled","Sleaf.modelled"))
error = as.factor(c("Mleaf_SD","Mstem_SD","Mroot_SD","Sleaf_SD"))
pd <- position_dodge(3) # move the overlapped errorbars horizontally
for (p in 1:length(meas)) {
  summary.data.Cpool = subset(summary.data,variable==meas[p])
  summary.output.Cpool = subset(summary.output,variable==res[p])
  summary.error.Cpool = subset(summary.error,variable==error[p])
  summary.error.Cpool$parameter = summary.data.Cpool$value
  p3 = ggplot(summary.error.Cpool, aes(x=Date, y=parameter, colour=volume, group=volume)) + 
    geom_errorbar(aes(ymin=parameter-value, ymax=parameter+value), colour="grey", width=3, position=pd) +
    geom_line(position=pd,data = summary.output.Cpool, aes(x = Date, y = value, group = interaction(volume,no.param), linetype=no.param, colour=volume)) + 
    geom_point(position=pd, size=2) +
    ylab(as.character(meas[p])) +
    ggtitle("C pools - Measured (points) vs Modelled (lines)") +
    labs(colour="Soil Volume", linetype="Total No of Parameter", shape="Total No of Parameter") +
    theme(legend.title = element_text(colour="chocolate", size=10, face="bold"))
  p3
  ggsave(p3,filename=paste(meas[p],"_Measured_vs_Modelled.png",sep=""))
}
plots5 <- lapply(ll <- list.files(patt='.*[.]png'),function(x){
  img <- as.raster(readPNG(x))
  rasterGrob(img, interpolate = FALSE)
})
ggsave("Summar_Cpools_multipage.pdf", marrangeGrob(grobs=plots5, nrow=2, ncol=2))


# Set working directory for saving figures
setwd("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/output/figures/summary/AF")

# Plot individual modelled parameters ("k","Y","af","as","ar","sf") against "volume" and "Total No of param"
var = as.factor(c("k","Y","af","as","ar","sf"))
for (p in 1:length(var)) {
  summary.param.set = subset(summary.param,variable==var[p])
  pd <- position_dodge(0.5) # move the overlapped errorbars horizontally
  p4 = ggplot() +
    geom_point(position=pd,data = summary.param.set, aes(x = Date, y = Parameter,  group = interaction(volume,no.param), colour=factor(volume), shape=factor(no.param))) +
    geom_line(position=pd,data = summary.param.set, aes(x = Date, y = Parameter,  group = interaction(volume,no.param), colour=factor(volume), linetype=factor(no.param))) +
    xlab("Days") +
    ylab(as.character(var[p])) +
    ggtitle("Modelled coefficients") +
    labs(colour="Soil Volume", linetype="Total No of Parameter", shape="Total No of Parameter") +
    theme(legend.title = element_text(colour="chocolate", size=10, face="bold"))
  ggsave(p4,filename=paste(var[p],"_over_time.png",sep=""))
  p4
}
plots6 <- lapply(ll <- list.files(patt='.*[.]png'),function(x){
  img <- as.raster(readPNG(x))
  rasterGrob(img, interpolate = FALSE)
})
ggsave("Summary_AF_multipage.pdf", marrangeGrob(grobs=plots6, nrow=2, ncol=2))


# Set working directory for saving figures
setwd("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/output/figures/summary")

# Plot modelled Cstorage against "volume" and "Total No of param"
p5 = ggplot() +
  geom_line(data = summary.Cstorage, aes(x = Date, y = Cstorage.modelled, group = interaction(volume,no.param),colour=volume, linetype=no.param)) + 
  ylab("Cstorage (gC)") +
  ggtitle("Modelled Cstorage") +
  labs(colour="Soil Volume", linetype="Total No of Parameter") +
  theme(legend.title = element_text(colour="chocolate", size=10, face="bold"))
p5
ggsave(p5,filename=paste("Cstorage_Modelled.png",sep=""))
ggsave(p5,filename=paste("Cstorage_Modelled.eps",sep=""))


# Plot modelled parameter means ("k","Y","af","as","ar","sf") against "volume" and "Total No of param"
param.data = as.data.frame(t(param.data))
names(param.data) = c("k","Y","af","as","ar","sf")
param.data$volume = as.factor("Synthetic Case")
param.data$no.param = as.factor(1)
param.data$type = as.factor("true")

param.mean$type = as.factor("modelled")
param.mean = rbind(param.mean,param.data)
melted.param.mean = melt(param.mean, id.vars=c("no.param","volume","type"))

pd <- position_dodge(0.3)
p6 = ggplot(data = melted.param.mean, aes(x = variable, y = value, group = interaction(volume,no.param,type), shape=factor(no.param), colour=factor(type), fill=factor(volume))) +
  geom_point(position=pd, size=4) +
  xlab("Allocation fractions") +
  ylab("Value of the coefficients") +
  ggtitle("Modelled mean allocation fractions") +
  labs(shape="Total No of Parameter", colour="Type", fill="Soil Volume") +
  theme(legend.title = element_text(colour="chocolate", size=10, face="bold"))
p6
ggsave(p6,filename=paste("Modelled_mean_allocation_fractions.png"))
ggsave(p6,filename=paste("Modelled_mean_allocation_fractions.eps"))


# Plot Model Measures ("logLi","aic","bic","time") against "volume" and "Total No of param"
pd <- position_dodge(0.3)
p7 = ggplot(data = melted.aic.bic, aes(x = variable, y = value, group = interaction(volume,no.param), shape=factor(no.param), colour=factor(volume))) +
  geom_point(position=pd, size=4) +
  xlab("Model Measures") +
  ylab("LogLi, AIC, BIC, Time") +
  ggtitle("LogLi, AIC, BIC, Time for various models") +
  labs(colour="Soil Volume", shape="Total No of Parameter") +
  theme(legend.title = element_text(colour="chocolate", size=10, face="bold"))
p7
ggsave(p7,filename=paste("LogLi_aic_bic_time.png"))
ggsave(p7,filename=paste("LogLi_aic_bic_time.eps"))

plots7 <- lapply(ll <- list.files(patt='.*[.]png'),function(x){
  img <- as.raster(readPNG(x))
  rasterGrob(img, interpolate = FALSE)
})
ggsave("Summary_rest_multipage.pdf", marrangeGrob(grobs=plots7, nrow=3, ncol=1))

# png("test.png")
# multiplot(p5, p6, p6, cols=1)
# dev.off()

