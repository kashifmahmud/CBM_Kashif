# Carbon balance model 
# Developed by Kashif Mahmud and Belinda Medlyn (November 2016)
# k.mahmud@westernsydney.edu.au

# This script creates the figures and saves those
##############################

# Set working directory for saving figures
setwd("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/output/figures")


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
  p3 = ggplot(summary.error.Cpool, aes(x=Date, y=parameter, shape=volume, group=volume)) + 
    geom_errorbar(aes(ymin=parameter-value, ymax=parameter+value), colour="grey", width=3, position=pd) +
    geom_line(position=pd,data = summary.output.Cpool, aes(x = Date, y = value, group = interaction(volume,no.param), linetype=volume, colour=no.param)) + 
    geom_point(position=pd, size=2) +
    ylab(as.character(meas[p])) +
    ggtitle("C pools - Measured (points) vs Modelled (lines)") +
    labs(shape="Soil Volume", linetype="Soil Volume", colour="Total No of Parameter") +
    theme(legend.title = element_text(colour="chocolate", size=10, face="bold"))
  p3
  ggsave(p3,filename=paste(meas[p],"_Measured_vs_Modelled.png",sep=""))
}


# Plot modelled Cstorage against "volume" and "Total No of param"
p4 = ggplot() +
  geom_line(data = summary.Cstorage, aes(x = Date, y = Cstorage.modelled, group = interaction(volume,no.param),colour=factor(no.param), linetype=factor(volume))) + 
  ylab("Cstorage (gC)") +
  ggtitle("Modelled Cstorage") +
  labs(linetype="Soil Volume", colour="Total No of Parameter") +
  theme(legend.title = element_text(colour="chocolate", size=10, face="bold"))
p4
ggsave(p4,filename=paste("Cstorage_Modelled.png",sep=""))


# Plot individual modelled parameters ("k","Y","af","as","ar","sf") against "volume" and "Total No of param"
var = as.factor(c("k","Y","af","as","ar","sf"))
for (p in 1:length(var)) {
  summary.param.set = subset(summary.param,variable==var[p])
  pd <- position_dodge(0.5) # move the overlapped errorbars horizontally
  p5 = ggplot() +
    geom_point(position=pd,data = summary.param.set, aes(x = Date, y = Parameter,  group = interaction(volume,no.param), colour=factor(no.param), shape=factor(volume))) +
    geom_line(position=pd,data = summary.param.set, aes(x = Date, y = Parameter,  group = interaction(volume,no.param), colour=factor(no.param), linetype=factor(volume))) +
    xlab("Days") +
    ylab(as.character(var[p])) +
    ggtitle("Modelled coefficients") +
    labs(shape="Soil Volume", linetype="Soil Volume", colour="Total No of Parameter") +
    theme(legend.title = element_text(colour="chocolate", size=10, face="bold"))
  ggsave(p5,filename=paste(var[p],"_over_time.png",sep=""))
  p5
}


# Plot modelled parameter means ("k","Y","af","as","ar","sf") against "volume" and "Total No of param"
melted.param.mean = melt(param.mean, id.vars=c("no.param","volume"))
pd <- position_dodge(0.1)
p6 = ggplot() +
  geom_point(position=pd, data = melted.param.mean, aes(x = variable, y = value, group = interaction(volume,no.param), colour=factor(no.param), shape=factor(volume))) +
  xlab("Allocation fractions") +
  ylab("Value of the coefficients") +
  ggtitle("Modelled allocation fractions") +
  labs(shape="Soil Volume", colour="Total No of Parameter") +
  theme(legend.title = element_text(colour="chocolate", size=10, face="bold"))
p6
ggsave(p6,filename=paste("Modelled_mean_allocation_fractions.png"))


# Plot Model Measures ("logLi","aic","bic","time") against "volume" and "Total No of param"
pd <- position_dodge(0.1)
p7 = ggplot(data = melted.aic.bic, aes(x = variable, y = value, group = interaction(volume,no.param), colour=factor(no.param), shape=factor(volume))) +
  geom_point(position=pd) +
  xlab("Model Measures") +
  ylab("LogLi, AIC, BIC, Time") +
  ggtitle("LogLi, AIC, BIC, Time for various models") +
  labs(shape="Soil Volume", colour="Total No of Parameter") +
  theme(legend.title = element_text(colour="chocolate", size=10, face="bold"))
p7
ggsave(p7,filename=paste("LogLi_aic_bic_time.png"))

