# Carbon balance model 
# Developed by Kashif Mahmud and Belinda Medlyn (November 2016)
# k.mahmud@westernsydney.edu.au

# This script creates the figures and saves those
##############################

# Set working directory for saving figures
setwd("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/output/figures/summary/AF")

# Plot individual modelled parameters ("k","Y","af","as","ar","sf") against "volume" and "Total No of param"
var = as.factor(c("k","Y","af","as","ar","sf"))
for (p in 1:length(var)) {
  summary.param.set.limit = subset(summary.param, variable==var[p])
  for (z in 1:length(no.param.par.var)) {
    summary.param.set = subset(summary.param, variable==var[p] & no.param==no.param.par.var[z])
    # summary.error.set = subset(summary.error, variable==var[p] & no.param==no.param.par.var[z])
    pd <- position_dodge(0.5) # move the overlapped errorbars horizontally
    p4 = ggplot(data = summary.param.set, aes(x = Date, y = Parameter,  group = volume.group, colour=factor(volume.group))) +
      geom_errorbar(data = summary.param.set, aes(ymin=Parameter-Parameter_SD, ymax=Parameter+Parameter_SD), width=0.5, size=0.1) +
      geom_point(position=pd,size=0.01) +
      geom_line(position=pd,data = summary.param.set, aes(x = Date, y = Parameter,  group = volume.group, colour=factor(volume.group))) +
      xlab("Days") +
      ylab(as.character(var[p])) +
      ggtitle(paste("Modelled coefficient,",as.character(var[p]),"for parameter",no.param.par.var[z])) +
      labs(colour="Treatment Group") +
      # scale_y_continuous(limits = c(param[1,1+(p-1)*3],param[1,3+(p-1)*3])) +
      scale_y_continuous(limits = c(min(summary.param.set.limit$Parameter)-max(summary.param.set.limit$Parameter_SD),max(summary.param.set.limit$Parameter)+max(summary.param.set.limit$Parameter_SD))) +
      # annotate("text", x = mean(summary.param.set.limit$Date), y = min(summary.param.set.limit$Parameter)-(mean(summary.param.set.limit$Parameter_SD)/4), size = 3,
      annotate("text", x = mean(summary.param.set.limit$Date), y = min(summary.param.set.limit$Parameter)+(mean(summary.param.set.limit$Parameter_SD)*2), size = 3,
                        label = paste("Group 1 = Volume: ", subset(summary.param.set.limit, volume.group==1)[1,5], "L", 
                             "\nGroup 2 = Volume: ", subset(summary.param.set.limit, volume.group==2)[1,5], "L",
                             "\nGroup 3 = Volume: ", subset(summary.param.set.limit, volume.group==3)[1,5], "L",
                             "\nGroup 4 = Volume: ", subset(summary.param.set.limit, volume.group==4)[1,5], "L",
                             # "\nGroup 5 = Volume: ", subset(summary.param.set.limit, volume.group==5)[1,5], "L",
                             # "\nGroup 6 = Volume: ", subset(summary.param.set.limit, volume.group==6)[1,5], "L",
                             # "\nGroup 7 = Volume: ", subset(summary.param.set.limit, volume.group==7)[1,5], "L",
                             # "\nGroup 8 = Volume: ", subset(summary.param.set.limit, volume.group==8)[1,5], "L",
                             "\nChain length = ", chainLength-bunr_in)) +
      theme_bw() +
      theme(plot.title = element_text(size = 12, face = "bold")) +
      theme(legend.title = element_text(colour="chocolate", size=12, face="bold")) +
      theme(axis.title.x = element_text(size = 12, vjust=-.2)) +
      theme(axis.title.y = element_text(size = 12, vjust=0.3))
    p4
    ggsave(p4,filename=paste(var[p],"_over_time_par_",z,".png",sep=""))
  }
}
plots6 <- lapply(ll <- list.files(patt='.*[.]png'),function(x){
  img <- as.raster(readPNG(x))
  rasterGrob(img, interpolate = FALSE)
})
ggsave("Summary_AF_multipage.pdf", marrangeGrob(grobs=plots6, nrow=2, ncol=length(no.param.par.var)))


# Set working directory for saving figures
setwd("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/output/figures/summary/Cpools")

# Plot modelled vs measured data ("Mleaf","Mstem","Mroot","Sleaf") against "volume" and "Total No of param"
meas = as.factor(c("Mleaf","Mstem","Mroot","Sleaf"))
res = as.factor(c("Mleaf.modelled","Mstem.modelled","Mroot.modelled","Sleaf.modelled"))
error = as.factor(c("Mleaf_SD","Mstem_SD","Mroot_SD","Sleaf_SD"))
pd <- position_dodge(2) # move the overlapped errorbars horizontally
for (p in 1:length(meas)) {
  summary.data.Cpool = subset(summary.data,variable==meas[p])
  summary.output.Cpool = subset(summary.output,variable==res[p])
  summary.error.Cpool = subset(summary.error,variable==error[p])
  # summary.error.Cpool$parameter = summary.data.Cpool$value
  
  p3 = ggplot(summary.error.Cpool, aes(x=Date, y=parameter, group = volume, colour=volume)) + 
    geom_point(position=pd) +
    geom_errorbar(position=pd,aes(ymin=parameter-value, ymax=parameter+value), colour="grey", width=2) +
    # geom_line(position=pd,data = summary.output.Cpool, aes(x = Date, y = value, group = interaction(volume,volume.group,no.param), linetype=volume.group, colour=volume, size=no.param)) + 
    geom_line(position=pd,data = summary.output.Cpool, aes(x = Date, y = value, group = interaction(volume,no.param), linetype=no.param, colour=volume)) + 
    ylab(paste(as.character(meas[p]),"(g C) in log scale")) +
    ggtitle("C pools - Measured (points) vs Modelled (lines)") +
    labs(colour="Soil Volume", linetype="No. of Parameters") +
    # labs(colour="Soil Volume", linetype="Grouping treatment", size="Total No of Parameter") +
    # scale_color_manual(labels = c("Individuals", "One Group"), values = c("blue", "red")) +
    theme_bw() +
    # annotate("text", x = mean(summary.param.set.limit$Date), y = min(summary.param.set.limit$Parameter)-mean(summary.param.set.limit$Parameter_SD), size = 3,
             # label = paste("Grouping treatment 1 = Individual parameter sets for different treatments",
             #               "\nGrouping treatment 2 = One single parameter set for all treatments",
             #               "\nChain length = ", chainLength-bunr_in)) +
    scale_y_log10() +
    theme(plot.title = element_text(size = 12, face = "bold")) +
    theme(legend.title = element_text(colour="chocolate", size=12, face="bold")) +
    theme(axis.title.x = element_text(size = 12, vjust=-.2)) +
    theme(axis.title.y = element_text(size = 12, vjust=0.3))
  p3
  ggsave(p3,filename=paste(meas[p],"_Measured_vs_Modelled.png",sep=""))
}
plots5 <- lapply(ll <- list.files(patt='.*[.]png'),function(x){
  img <- as.raster(readPNG(x))
  rasterGrob(img, interpolate = FALSE)
})
ggsave("Summary_Cpools_multipage.pdf", marrangeGrob(grobs=plots5, nrow=2, ncol=2))


# Set working directory for saving figures
setwd("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/output/figures/summary")

# Plot modelled Cstorage against "volume" and "Total No of param"
p5 = ggplot() +
  # geom_line(data = summary.Cstorage, aes(x = Date, y = Cstorage.modelled, group = interaction(volume,volume.group),colour=volume, linetype=volume.group)) + 
  geom_line(data = summary.Cstorage, aes(x = Date, y = Cstorage.modelled, group = interaction(volume,no.param),colour=volume, linetype=no.param)) + 
  ylab("Cstorage (gC)") +
  ggtitle("Modelled Cstorage") +
  labs(colour="Soil Volume", linetype="Total No of Parameter") +
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold")) +
  theme(legend.title = element_text(colour="chocolate", size=12, face="bold")) +
  theme(axis.title.x = element_text(size = 12, vjust=-.2)) +
  theme(axis.title.y = element_text(size = 12, vjust=0.3))
p5
ggsave(p5,filename=paste("Cstorage_Modelled.png",sep=""))


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
p7 = ggplot(data = melted.aic.bic, aes(x = variable, y = value, group = interaction(volume.group,no.param), shape=factor(no.param), colour=factor(volume.group))) +
  geom_point(position=pd, size=4) +
  xlab("Model Measures") +
  ylab("LogLi, AIC, BIC, Time") +
  ggtitle("LogLi, AIC, BIC, Time for various models") +
  labs(colour="Treatment group", shape="Total No of Parameter") +
  annotate("text", x = melted.aic.bic$variable[1+nrow(melted.aic.bic)/2], y = min(melted.aic.bic$value)/2, size = 3,
           label = paste("Group 1 = Volume: ", subset(summary.param.set.limit, volume.group==1)[1,5], "L", 
                         "\nGroup 2 = Volume: ", subset(summary.param.set.limit, volume.group==2)[1,5], "L",
                         "\nGroup 3 = Volume: ", subset(summary.param.set.limit, volume.group==3)[1,5], "L",
                         "\nGroup 4 = Volume: ", subset(summary.param.set.limit, volume.group==4)[1,5], "L",
                         # "\nGroup 5 = Volume: ", subset(summary.param.set.limit, volume.group==5)[1,5], "L",
                         # "\nGroup 6 = Volume: ", subset(summary.param.set.limit, volume.group==6)[1,5], "L",
                         # "\nGroup 7 = Volume: ", subset(summary.param.set.limit, volume.group==7)[1,5], "L",
                         # "\nGroup 8 = Volume: ", subset(summary.param.set.limit, volume.group==8)[1,5], "L",
                         "\nChain length = ", chainLength-bunr_in)) +
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold")) +
  theme(legend.title = element_text(colour="chocolate", size=12, face="bold")) +
  theme(axis.title.x = element_text(size = 12, vjust=-.2)) +
  theme(axis.title.y = element_text(size = 12, vjust=0.3))
p7
ggsave(p7,filename=paste("LogLi_aic_bic_time.png"))

# # Plot Model Measures ("logLi","aic","bic","time") against "volume" and "Total No of param"
# pd <- position_dodge(0)
# p8 = ggplot(data = aic.bic[c("bic","no.param","volume")], aes(x = no.param, y = bic, group = volume, colour=factor(volume))) +
#   geom_line() +
#   geom_point(position=pd, size=2) +
#   xlab("Model Measures") +
#   ylab("BIC") +
#   ggtitle("BIC for various models") +
#   labs(colour="No. of Parameter") +
#   annotate("text", x = mean(melted.aic.bic$no.param), y = min(aic.bic$bic), size = 3,
#            label = paste("Group 1 = Volume: ", subset(summary.param.set.limit, volume.group==1)[1,5], "L", 
#                          "\nGroup 2 = Volume: ", subset(summary.param.set.limit, volume.group==2)[1,5], "L",
#                          # "\nGroup 3 = Volume: ", subset(summary.param.set.limit, volume.group==3)[1,5], "L",
#                          # "\nGroup 4 = Volume: ", subset(summary.param.set.limit, volume.group==4)[1,5], "L",
#                          # "\nGroup 5 = Volume: ", subset(summary.param.set.limit, volume.group==5)[1,5], "L",
#                          # "\nGroup 6 = Volume: ", subset(summary.param.set.limit, volume.group==6)[1,5], "L",
#                          # "\nGroup 7 = Volume: ", subset(summary.param.set.limit, volume.group==7)[1,5], "L",
#                          # "\nGroup 8 = Volume: ", subset(summary.param.set.limit, volume.group==8)[1,5], "L",
#                          "\nChain length = ", chainLength-bunr_in)) +
#   theme_bw() +
#   theme(plot.title = element_text(size = 12, face = "bold")) +
#   theme(legend.title = element_text(colour="chocolate", size=12, face="bold")) +
#   theme(axis.title.x = element_text(size = 12, vjust=-.2)) +
#   theme(axis.title.y = element_text(size = 12, vjust=0.3))
# p8
# ggsave(p8,filename=paste("bic_with_storage.png"))

# # Plot Model Measures ("bic") against "volume.group" and "volume"
# bic.final = read.csv("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/output/processeddata/bic.final.csv")
# bic.final$volume = as.factor(bic.final$volume)
# pd <- position_dodge(0.4)
# p9 = ggplot(data = bic.final, aes(x = volume, y = bic, group = volume.group, colour=factor(volume.group))) +
#   # geom_line() +
#   geom_point(position=pd, size=2) +
#   xlab("Pot volume (L)") +
#   ylab("BIC") +
#   ggtitle("BIC for various parameter settings") +
#   labs(colour="Volume group") +
#   annotate("text", x = bic.final$volume[4], y = max(bic.final$bic)-400, size = 3,
#            label = paste("Group 1 = c(5,10,15,20,25,35,1000L)",
#                          "\nGroup 2 = c(5,10,15,20,25,35L) and 1000L",
#                          "\nGroup 3 = c(5,10L), c(15,20,25,35L) and 1000L",
#                          "\nGroup 4 = c(5,10L), c(15,20L), c(25,35L) and 1000L",
#                          "\nGroup 5 = Individual treatment parameters",
#                          # "\nGroup 6 = Volume: ", subset(summary.param.set.limit, volume.group==6)[1,5], "L",
#                          # "\nGroup 7 = Volume: ", subset(summary.param.set.limit, volume.group==7)[1,5], "L",
#                          # "\nGroup 8 = Volume: ", subset(summary.param.set.limit, volume.group==8)[1,5], "L",
#                          "\nChain length = ", chainLength-bunr_in)) +
#   theme_bw() +
#   theme(plot.title = element_text(size = 12, face = "bold")) +
#   theme(legend.title = element_text(colour="chocolate", size=12, face="bold")) +
#   theme(axis.title.x = element_text(size = 12, vjust=-.2)) +
#   theme(axis.title.y = element_text(size = 12, vjust=0.3))
# p9
# # ggsave(p9,filename=paste("bic_with_storage.png"))

plots7 <- lapply(ll <- list.files(patt='.*[.]png'),function(x){
  img <- as.raster(readPNG(x))
  rasterGrob(img, interpolate = FALSE)
})
ggsave("Summary_rest_multipage.pdf", marrangeGrob(grobs=plots7, nrow=2, ncol=1))

# png("test.png")
# multiplot(p5, p6, p6, cols=1)
# dev.off()

