# Figure plotting for Paper 01 (CDM with MCMC for sink limited pot experiment)

# Set working directory for saving figures
setwd("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/archive/processeddata")

################# Figure 3 ####################
# Plot Model Measures ("bic") against "models with and without storage pool"
bic.with.storage = read.csv("logli_aic_bic_time_with_storage.csv")
bic.without.storage = read.csv("logli_aic_bic_time_without_storage.csv")
keeps = c("bic", "volume")
bic.with.storage = bic.with.storage[ , keeps, drop = FALSE]
bic.without.storage = bic.without.storage[ , keeps, drop = FALSE]
bic = merge(bic.with.storage,bic.without.storage,by=c("volume"))
names(bic)[2:3] <- c("bic with storage", "bic without storage")
bic$volume = as.factor(bic$volume)
bic.melt <- melt(bic, id.vars = "volume")

pd <- position_dodge(0)
p1 = ggplot(data = bic.melt, aes(x = volume, y = value, group = variable, shape=factor(variable))) +
  geom_line(data = bic.melt, aes(x = volume, y = value, group = variable, colour=factor(variable)), show.legend=FALSE) +
  geom_point(position=pd, size=3) +
  xlab("Pot volume (L)") +
  ylab("BIC") +
  ggtitle("BIC for various model settings") +
  labs(shape="Model setting") +
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold")) +
  theme(legend.title = element_text(colour="chocolate", size=12, face="bold")) +
  theme(legend.position = c(0.2,0.8)) +
  theme(legend.key = element_blank()) +
  theme(axis.title.x = element_text(size = 12, vjust=-.2)) +
  theme(axis.title.y = element_text(size = 12, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p1
ggsave(p1,filename=paste("bic_with_without_storage.png"))


################ Figure 4 #####################
# Plot Model Measures ("bic") against "Treatment groupings"
bic.final = read.csv("bic.all.csv")
# bic.final = read.csv("bic.final.csv")
bic.final$volume = as.factor(bic.final$volume)

pd <- position_dodge(0.3)
p2 = ggplot(data = bic.final, aes(x = volume, y = bic, group = volume.group, colour=factor(volume.group))) +
  geom_line(position=pd, show.legend=FALSE) +
  geom_point(position=pd, size=2) +
  xlab("Pot volume (L)") +
  ylab("BIC") +
  ggtitle("BIC for various treatment grouping") +
  labs(colour="Volume group") +
  # annotate("text", x = bic.final$volume[4], y = max(bic.final$bic)-400, size = 3,
  #          label = paste("Group 1 = c(5,10,15,20,25,35,1000L)",
  #                        "\nGroup 2 = c(5,10,15,20,25,35L) and 1000L",
  #                        "\nGroup 3 = c(5,10L), c(15,20,25,35L) and 1000L",
  #                        "\nGroup 4 = c(5,10L), c(15,20L), c(25,35L) and 1000L",
  #                        "\nGroup 5 = Individual treatment parameters",
  #                        # "\nGroup 6 = Volume: ", subset(summary.param.set.limit, volume.group==6)[1,5], "L",
  #                        # "\nGroup 7 = Volume: ", subset(summary.param.set.limit, volume.group==7)[1,5], "L",
  #                        # "\nGroup 8 = Volume: ", subset(summary.param.set.limit, volume.group==8)[1,5], "L",
  #                        "\nChain length = ", chainLength-bunr_in)) +
  scale_colour_discrete(name="Treatment group",
                      breaks=c("1", "2", "3", "4", "5", "6", "7"),
                      labels=c("Group 1 = c(5,10,15,20,25,35,1000L)", "Group 2 = c(5,10,15,20,25,35L) and 1000L", 
                               "Group 3 = c(5,10L), c(15,20,25,35L) and 1000L", "Group 4 = c(5,10L), c(15,20L), c(25,35L) and 1000L",
                               "Group 5 = c(5,10L), c(15,20,25), 35L and 1000L", "Group 6 = c(5,10L), 20L, c(15,25,35L) and 1000L",
                               "Group 7 = Individual treatment parameters")) +
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold")) +
  # theme(legend.title = element_text(colour="chocolate", size=12, face="bold")) +
  theme(legend.text = element_text(colour="blue", size = 8)) +
  theme(legend.key.height=unit(0.7,"line")) +
  theme(legend.position = c(0.3,0.78)) +
  theme(legend.key = element_blank()) +
  theme(axis.title.x = element_text(size = 12, vjust=-.2)) +
  theme(axis.title.y = element_text(size = 12, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p2
ggsave(p2,filename=paste("bic_treat_group.png"))


################ Figure 5 #####################
# Plot Model Measures ("bic") against "Total No of parameter"
bic.param = read.csv("logli_aic_bic_time.csv")
keeps = c("bic", "volume.group", "no.param")
bic.param = bic.param[ , keeps, drop = FALSE]
  
pd <- position_dodge(0)
p3 = ggplot(data = bic.param[c("bic","no.param","volume.group")], aes(x = no.param, y = bic, group = volume.group, colour=factor(volume.group))) +
  geom_line() +
  geom_point(position=pd, size=2) +
  xlab("Number of parameter") +
  ylab("BIC") +
  ggtitle("BIC for various parameter numbers") +
  labs(colour="Treatment groups") +
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold")) +
  theme(legend.title = element_text(colour="chocolate", size=12, face="bold")) +
  theme(axis.title.x = element_text(size = 12, vjust=-.2)) +
  theme(axis.title.y = element_text(size = 12, vjust=0.3))
p3
ggsave(p3,filename=paste("bic_param_number.png"))

