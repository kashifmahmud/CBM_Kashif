# Carbon balance model 
# Developed by Kashif Mahmud and Belinda Medlyn (November 2016)
# k.mahmud@westernsydney.edu.au

# This script creates the figures and saves those
##############################
# This version is for previous model without storage pool

# Set working directory for saving figures
setwd("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/output/figures/summary/Cpools")

# Plot modelled vs measured data ("Mleaf","Mstem","Mroot","Sleaf") against "volume" and "Total No of param"
meas = as.factor(c("Mleaf","Mstem","Mroot"))
res = as.factor(c("Mleaf.modelled","Mstem.modelled","Mroot.modelled"))
error = as.factor(c("Mleaf_SD","Mstem_SD","Mroot_SD"))
pd <- position_dodge(2) # move the overlapped errorbars horizontally
for (p in 1:length(meas)) {
  summary.data.Cpool = subset(summary.data,variable==meas[p])
  summary.output.Cpool = subset(summary.output,variable==res[p])
  summary.error.Cpool = subset(summary.error,variable==error[p])
  # summary.error.Cpool$parameter = summary.data.Cpool$value
  
  p3 = ggplot(summary.error.Cpool, aes(x=Date, y=parameter, group = volume, colour=volume)) + 
    geom_point(position=pd) +
    geom_errorbar(aes(ymin=parameter-value, ymax=parameter+value), colour="grey", width=2) +
    geom_line(position=pd,data = summary.output.Cpool, aes(x = Date, y = value, group = interaction(volume,no.param), linetype=no.param, colour=volume)) + 
    ylab(paste(as.character(meas[p]),"(g C) in log scale")) +
    ggtitle("C pools - Measured (points) vs Modelled (lines) without storage") +
    labs(colour="Soil Volume", linetype="Total No of Parameter", shape="Total No of Parameter") +
    scale_y_log10() +
    theme_bw() +
    theme(plot.title = element_text(size = 12, face = "bold")) +
    theme(legend.title = element_text(colour="chocolate", size=12, face="bold")) +
    theme(axis.title.x = element_text(size = 12, vjust=-.2)) +
    theme(axis.title.y = element_text(size = 12, vjust=0.3))
  p3
  ggsave(p3,filename=paste(meas[p],"_Measured_vs_Modelled_without_storage.png",sep=""))
}
plots5 <- lapply(ll <- list.files(patt='.*[.]png'),function(x){
  img <- as.raster(readPNG(x))
  rasterGrob(img, interpolate = FALSE)
})
ggsave("Summary_Cpools_multipage.pdf", marrangeGrob(grobs=plots5, nrow=2, ncol=2))


# Set working directory for saving figures
setwd("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/output/figures/summary/AF")

# Plot individual modelled parameters ("k","Y","af","as","ar","sf") against "volume" and "Total No of param"
var = as.factor(c("Y","af","as","ar","sf"))
for (p in 1:length(var)) {
  summary.param.set.limit = subset(summary.param, variable==var[p])
  for (z in 1:length(no.param.par.var)) {
    summary.param.set = subset(summary.param, variable==var[p] & no.param==no.param.par.var[z])
    # summary.error.set = subset(summary.error, variable==var[p] & no.param==no.param.par.var[z])
    pd <- position_dodge(0.5) # move the overlapped errorbars horizontally
    p4 = ggplot(data = summary.param.set, aes(x = Date, y = Parameter,  group = volume, colour=factor(volume))) +
      geom_errorbar(data = summary.param.set, aes(ymin=Parameter-Parameter_SD, ymax=Parameter+Parameter_SD), width=0.5, size=0.1) +
      geom_point(position=pd,size=0.01) +
      geom_line(position=pd,data = summary.param.set, aes(x = Date, y = Parameter,  group = volume, colour=factor(volume))) +
      xlab("Days") +
      ylab(as.character(var[p])) +
      ggtitle(paste("Modelled coefficient,",as.character(var[p]),"for parameter",no.param.par.var[z],"without storage")) +
      labs(colour="Soil Volume") +
      # scale_y_continuous(limits = c(param[1,1+(p-1)*3],param[1,3+(p-1)*3])) +
      # scale_y_continuous(limits = c(min(summary.param.set.limit$Parameter)-0.25,max(summary.param.set.limit$Parameter)+0.25)) +
      scale_y_continuous(limits = c(min(summary.param.set.limit$Parameter)-max(summary.param.set.limit$Parameter_SD),max(summary.param.set.limit$Parameter)+max(summary.param.set.limit$Parameter_SD))) +
      theme_bw() +
      theme(plot.title = element_text(size = 12, face = "bold")) +
      theme(legend.title = element_text(colour="chocolate", size=12, face="bold")) +
      theme(axis.title.x = element_text(size = 12, vjust=-.2)) +
      theme(axis.title.y = element_text(size = 12, vjust=0.3))
    p4
    ggsave(p4,filename=paste(var[p],"_over_time_par_",z,"_without_storage",".png",sep=""))
  }
}
plots6 <- lapply(ll <- list.files(patt='.*[.]png'),function(x){
  img <- as.raster(readPNG(x))
  rasterGrob(img, interpolate = FALSE)
})
# ggsave("Summary_AF_multipage.pdf", marrangeGrob(grobs=plots6, nrow=2, ncol=length(no.param.par.var)))
ggsave("Summary_AF_multipage.pdf", marrangeGrob(grobs=plots6, nrow=2, ncol=2))


# Set working directory for saving figures
setwd("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/output/figures/summary")

# # Plot modelled Cstorage against "volume" and "Total No of param"
# p5 = ggplot() +
#   geom_line(data = summary.Cstorage, aes(x = Date, y = Cstorage.modelled, group = interaction(volume,no.param),colour=volume, linetype=no.param)) + 
#   ylab("Cstorage (gC)") +
#   ggtitle("Modelled Cstorage") +
#   labs(colour="Soil Volume", linetype="Total No of Parameter") +
#   theme_bw() +
#   theme(plot.title = element_text(size = 12, face = "bold")) +
#   theme(legend.title = element_text(colour="chocolate", size=12, face="bold")) +
#   theme(axis.title.x = element_text(size = 12, vjust=-.2)) +
#   theme(axis.title.y = element_text(size = 12, vjust=0.3))
# p5
# ggsave(p5,filename=paste("Cstorage_Modelled.png",sep=""))


# # Plot modelled parameter means ("k","Y","af","as","ar","sf") against "volume" and "Total No of param"
# melted.param.mean = melt(param.mean, id.vars=c("no.param","volume"))
# pd <- position_dodge(0.3)
# p6 = ggplot(data = melted.param.mean, aes(x = variable, y = value, group = interaction(volume,no.param), shape=factor(no.param), colour=factor(volume))) +
#   geom_point(position=pd, size=4) +
#   xlab("Allocation fractions") +
#   ylab("Value of the coefficients") +
#   ggtitle("Modelled mean allocation fractions") +
#   labs(colour="Soil Volume", shape="Total No of Parameter") +
#   theme_bw() +
#   theme(plot.title = element_text(size = 12, face = "bold")) +
#   theme(legend.title = element_text(colour="chocolate", size=12, face="bold")) +
#   theme(axis.title.x = element_text(size = 10, vjust=-.2)) +
#   theme(axis.title.y = element_text(size = 10, vjust=0.3))
# p6
# ggsave(p6,filename=paste("Modelled_mean_allocation_fractions.png"))


# Plot Model Measures ("logLi","aic","bic","time") against "volume" and "Total No of param"
pd <- position_dodge(0.3)
p7 = ggplot(data = melted.aic.bic, aes(x = variable, y = value, group = interaction(volume,no.param), shape=factor(no.param), colour=factor(volume))) +
  geom_point(position=pd, size=4) +
  xlab("Model Measures") +
  ylab("LogLi, AIC, BIC, Time") +
  ggtitle("LogLi, AIC, BIC, Time for various models without storage") +
  labs(colour="Soil Volume", shape="Total No of Parameter") +
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold")) +
  theme(legend.title = element_text(colour="chocolate", size=12, face="bold")) +
  theme(axis.title.x = element_text(size = 12, vjust=-.2)) +
  theme(axis.title.y = element_text(size = 12, vjust=0.3))
p7
ggsave(p7,filename=paste("LogLi_aic_bic_time_without_storage.png"))

plots7 <- lapply(ll <- list.files(patt='.*[.]png'),function(x){
  img <- as.raster(readPNG(x))
  rasterGrob(img, interpolate = FALSE)
})
ggsave("Summary_rest_multipage.pdf", marrangeGrob(grobs=plots7, nrow=1, ncol=2))

# # Plot Model Measures ("bic") against "models" and "treatments"
# bic.with.storage = read.csv("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/output/processeddata/logli_aic_bic_time_with_storage.csv")
# bic.without.storage = read.csv("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/output/processeddata/logli_aic_bic_time_without_storage.csv")
# keeps <- c("bic", "volume")
# bic.with.storage = bic.with.storage[ , keeps, drop = FALSE]
# bic.without.storage = bic.without.storage[ , keeps, drop = FALSE]
# bic = merge(bic.with.storage,bic.without.storage,by=c("volume"))
# names(bic)[2:3] <- c("bic with storage", "bic without storage")
# bic$volume = as.factor(bic$volume)
# bic.melt <- melt(bic, id.vars = "volume")
# 
# pd <- position_dodge(0)
# p9 = ggplot(data = bic.melt, aes(x = volume, y = value, group = variable, shape=factor(variable))) +
#   geom_line(data = bic.melt, aes(x = volume, y = value, group = variable, colour=factor(variable)), show.legend=FALSE) +
#   geom_point(position=pd, size=3) +
#   xlab("Pot volume (L)") +
#   ylab("BIC") +
#   ggtitle("BIC for various model settings") +
#   labs(shape="Model setting") +
#   theme_bw() +
#   theme(plot.title = element_text(size = 12, face = "bold")) +
#   theme(legend.title = element_text(colour="chocolate", size=12, face="bold")) +
#   theme(legend.position = c(0.2,0.8)) +
#   theme(legend.key = element_blank()) +
#   theme(axis.title.x = element_text(size = 12, vjust=-.2)) +
#   theme(axis.title.y = element_text(size = 12, vjust=0.3))
# p9
# ggsave(p9,filename=paste("bic_with_without_storage.png"))


