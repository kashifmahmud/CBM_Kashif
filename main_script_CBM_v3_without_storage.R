# Carbon balance model (CBM)
# Developed by Kashif Mahmud and Belinda Medlyn (November 2016)
# k.mahmud@westernsydney.edu.au

# Updated version without treatment groupings

# This code carries out Bayesian calibration for 5 variables (allocation fractions: "k","Y",af","as","sf") on 
# various temporal scales (e.g. 1,2,...,121 days) to estimate Carbon pools (Cstorage,Cleaf,Cstem,Croot) and fluxes

##############################
# MCMC with soil manipulation pot experiment data for all treatments (including the free seedling), 
# This version considers either daily/weekly/monthly/just one parameter set for 5 variables ("k","Y","af","as","sf")
# So we can set the parameters for various time frames
# Also calculates the MCMC SDs for all parameters with different time frames, and also the LogLi, AIC, BIC, time (measures 
# for best model selection) to select the best parameter set
# Finally save the figures in Github/results folder
##############################
# Set seed number to reproduce the results
set.seed(1)

# Set working directory for saving figures
setwd("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif")

# This script cleans the workspace, loads necessary Rfunctions and packages
source("load_packages_functions_CBM.R")

# Load the function to define the CBM equations to iteratively calculate Cstorage, Cleaf, Cstem, Croot, Sleaf, Sstem, Sroot
source("Rfunctions/CBM_model_v3_without_storage.R")

# # This script imports and processes the raw HIE pot experiment data to model the carbon pools and fluxes using MCMC
# source("initial_data_processing.R")

# This sript reads the Pot experiment raw data
source("read_data_CBM.R")

# Assign inputs for MCMC
chainLength = 3000 # Setting the length of the Markov Chain to be generated
bunr_in = chainLength * 0.1 # Discard the first 10% iterations for Burn-IN in MCMC (According to Oijen, 2008)
no.var = 4 # variables to be modelled are: Y,af,as,sf (without storage pool, so no "k")

# Assign pot volumes and number of parameters per varible in temporal scale
# vol = c(20) # test run
no.param.par.var = c(1) # test run
GPP.data.raw = read.csv("rawdata/GPP.csv") # Units gC d-1
vol = unique(GPP.data.raw$volume)[order(unique(GPP.data.raw$volume))] # Assign all treatment pot volumes
# no.param.par.var = c(1,2,3,4,5,6,9) # temporal parameter count per variable

param.mean = data.frame(matrix(ncol = no.var+1, nrow = length(no.param.par.var)*length(vol)))
names(param.mean) = c("Y","af","as","ar","sf")
aic.bic = data.frame(matrix(ncol = 4, nrow = length(no.param.par.var)*length(vol)))
time = data.frame(no.param=rep(no.param.par.var,length(vol)),
                  start.time=numeric(length(no.param.par.var)*length(vol)),
                  end.time=numeric(length(no.param.par.var)*length(vol)),
                  time.taken=numeric(length(no.param.par.var)*length(vol)))
q = 0 # Indicates the iteration number

# z=1; v=1
for (v in 1:length(vol)) {
  for (z in 1:length(no.param.par.var)) {
    # for (v in 1:length(vol)) {
    # This script process the raw data
    source("data_processing_CBM_v3.R")
    
    
    # Initialize few output data files
    q = q + 1
    time$start.time[q] <- Sys.time()
    param.vary = ceiling(nrow(data)/no.param.par.var[z]) # How many days the parameter set remain unchanged (weekly = 7; monthly = 30; just one parameter = nrow(data))
    no.param = ceiling(nrow(data)/param.vary) # number of parameter set for the whole duration of experiment (121 days)
    # j = c()
    # j[1] = 0
    # i = seq(1,nrow(data),1)
    # j[i] = i - ceiling(i/param.vary)*1  # j is for parameter settings for various time frames
    # 
    
    # This script initializes the parameter setting
    source("parameter_setting_v3.R")
    
    # Defining the variance-covariance matrix for proposal generation
    vcov = (0.01*(pMaxima-pMinima))^2
    vcovProposal =  vcov # The higher the coefficient, the higher the deviations in parameter time series
    # vcovProposal =  vcov[1,] # The higher the coefficient, the higher the deviations in parameter time series
    # with lower acceptance rate and better matching
    
    
    # Find the Prior probability density
    prior.dist = vector("list", no.var)
    for (i in 1:no.var) {
      prior.dist[i] = list(log(dunif(pValues[ , i], pMinima[ , i], pMaxima[ , i])))
    }
    logPrior0 <- sum(unlist(prior.dist))
    
    
    # Calculating model outputs for the starting point of the chain
    Mleaf = Mstem = Mroot = c()
    Mleaf[1] <- data$Mleaf[1]
    Mstem[1] <- data$Mstem[1]
    Mroot[1] <- data$Mroot[1]
    
    # GPP=data$GPP; Rd=data$Rd
    # Y=pValues$Y; k=pValues$k; af=pValues$af; as=pValues$as; sf=pValues$sf
    
    output = model(data$GPP,data$Rd,no.param,Mleaf,Mstem,Mroot,pValues$Y,pValues$af,pValues$as,pValues$sf)
    
    # Modification to consider the mean sf values over all treatments 
    # output = model(data$GPP,data$Rd,no.param,Mleaf,Mstem,Mroot,pValues$Y,pValues$k,pValues$af,pValues$as,param.sf.mean$sf.mean)
    
    
    logL0 <- logLikelihood(data,output) # Calculate log likelihood of starting point of the chain
    pChain[1,] <- c(pValues$Y,pValues$af,pValues$as,pValues$sf,logL0) # Assign the first parameter set with log likelihood
    # pChain[1,] <- c(pValues$k,pValues$Y,pValues$af,pValues$as,param.sf.mean$sf.mean,logL0) # Assign the first parameter set with log likelihood
    
    
    # Calculating the next candidate parameter vector, as a multivariate normal jump away from the current point
    # c=2
    for (c in (2 : chainLength)) {
      candidatepValues = matrix(ncol = no.var, nrow = no.param)
      for (i in 1:no.var) {
        # for (j in 1:no.param) {
        candidatepValues[,i] = rmvnorm(n=1, mean=pValues[,i],
                                       sigma=diag(vcovProposal[,i],no.param)) 
        # }
      }
      candidatepValues = data.frame(candidatepValues)
      names(candidatepValues) <- c("Y","af","as","sf")
      
      
      # Reflected back to generate another candidate value
      reflectionFromMin = pmin( unlist(matrix(0,nrow=no.param,ncol=no.var)), unlist(candidatepValues-pMinima) )
      reflectionFromMax = pmax( unlist(list(rep(0, no.param))), unlist(candidatepValues-pMaxima) )
      candidatepValues = candidatepValues - 2 * reflectionFromMin - 2 * reflectionFromMax 
      
      
      # Calculating the prior probability density for the candidate parameter vector
      if (all(candidatepValues>pMinima) && all(candidatepValues<pMaxima)){
        uni.dist = vector("list", no.var)
        for (i in 1:no.var) {
          uni.dist[i] = list(log(dunif(candidatepValues[ , i], pMinima[ , i], pMaxima[ , i])))
        }
        logPrior1 <- sum(unlist(uni.dist))
        Prior1 = 1
      } else {
        Prior1 <- 0
      }
      
      
      # Calculating the outputs for the candidate parameter vector and then log likelihood
      if (Prior1 > 0) {
        Mleaf = Mstem = Mroot = c()
        Mleaf[1] <- data$Mleaf[1]
        Mstem[1] <- data$Mstem[1]
        Mroot[1] <- data$Mroot[1]
        
        GPP=data$GPP; Rd=data$Rd
        Y=candidatepValues$Y; af=candidatepValues$af; as=candidatepValues$as; sf=candidatepValues$sf
        out.cand = model(data$GPP,data$Rd,no.param,Mleaf,Mstem,Mroot,candidatepValues$Y,
                         candidatepValues$af,candidatepValues$as,candidatepValues$sf)
        
        # Modification to consider the mean sf values over all treatments
        # sf=param.sf.mean$sf.mean
        # out.cand = model(data$GPP,data$Rd,no.param,Mleaf,Mstem,Mroot,candidatepValues$Y,
        #                  candidatepValues$k,candidatepValues$af,candidatepValues$as,param.sf.mean$sf.mean)
        
        logL1 <- logLikelihood(data,out.cand) # Calculate log likelihood
        
        
        # Calculating the logarithm of the Metropolis ratio
        logalpha <- (logPrior1+logL1) - (logPrior0+logL0) 
        # Accepting or rejecting the candidate vector
        if ( log(runif(1, min = 0, max =1)) < logalpha ) {
        # if ( log(runif(1, min = 0, max =1)) < logalpha && candidatepValues$af + candidatepValues$as <= 1 
        #      && candidatepValues$as >= 0 && candidatepValues$af >= 0) {
          # && param.k[3] <= candidatepValues$k <= param.k[1,1] && param.Y[3] <= candidatepValues$Y <= param.Y[1] ) {
          pValues <- candidatepValues
          logPrior0 <- logPrior1
          logL0 <- logL1
        }
      }
      pChain[c,] <- c(pValues$Y,pValues$af,pValues$as,pValues$sf,logL0)
      
      # Modification to consider the mean sf values over all treatments 
      # pChain[c,] <- c(pValues$k,pValues$Y,pValues$af,pValues$as,param.sf.mean$sf.mean,logL0)
    }
    # Discard the first 500 iterations for Burn-IN in MCMC
    pChain <- pChain[(bunr_in+1):nrow(pChain),]
    pChain = as.data.frame(pChain)
    if (no.param.par.var[z]==1) {
      names(pChain) <- c("Y1","af1","as1","sf1","logli")
    }
    if (no.param.par.var[z]==2) {
      names(pChain) <- c("Y1","Y2","af1","af2","as1","as2","sf1","sf2","logli")
    }
    if (no.param.par.var[z]==3) {
      names(pChain) <- c("Y1","Y2","Y3","af1","af2","af3","as1","as2","as3","sf1","sf2","sf3","logli")
    }
    if (no.param.par.var[z]==4) {
      names(pChain) <- c("Y1","Y2","Y3","Y4","af1","af2","af3","af4","as1","as2","as3","as4","sf1","sf2","sf3","sf4","logli")
    }
    
    # Find the correlation matrix between parameter set
    corrMatrix = cor(pChain[,c(1:ncol(pChain)-1)])
    if (no.param.par.var[z]==1) {
      names(corrMatrix) <- c("Y1","af1","as1","sf1")
    }
    if (no.param.par.var[z]==2) {
      names(corrMatrix) <- c("Y1","Y2","af1","af2","as1","as2","sf1","sf2")
    }
    if (no.param.par.var[z]==3) {
      names(corrMatrix) <- c("Y1","Y2","Y3","af1","af2","af3","as1","as2","as3","sf1","sf2","sf3")
    }
    if (no.param.par.var[z]==4) {
      names(corrMatrix) <- c("Y1","Y2","Y3","Y4","af1","af2","af3","af4","as1","as2","as3","as4","sf1","sf2","sf3","sf4")
    }
    
    # Plotting the correlation matrix between parameter set
    png(height=1200, width=1200, pointsize=25, file = paste("output/figures/corrMatrix/corrMatrix_",v,"_vol_",vol[v],"_par_",no.param.par.var[z],".png",sep=""))
    corrplot(corrMatrix,tl.cex=1.5,title=paste("Correlation Matrix for vol",vol[v],"with par",no.param.par.var[z]), method="circle", is.corr=FALSE,type="full", cl.cex=2,
             addgrid.col="blue",addshade="positive", addCoef.col = rgb(0,0,0), mar=c(0,0,1,0), diag= FALSE,cl.lim = c(-1,1))
    dev.off()
    
    
    # Store the final parameter set values
    param.set = colMeans(pChain[ , 1:(no.param*no.var)])
    param.SD = apply(pChain[ , 1:(no.param*no.var)], 2, sd)
    param.final = data.frame(matrix(ncol = (no.var)*2, nrow = no.param))
    names(param.final) <- c("Y","af","as","sf","Y_SD","af_SD","as_SD","sf_SD")
    param.final$Y = param.set[1:no.param]
    param.final$af = param.set[(1+no.param):(2*no.param)]
    param.final$as = param.set[(1+2*no.param):(3*no.param)]
    param.final$sf = param.set[(1+3*no.param):(4*no.param)]
    # param.final$ar = 1 - param.final$af - param.final$as
    param.final$Y_SD = param.SD[1:no.param]
    param.final$af_SD = param.SD[(1+no.param):(2*no.param)]
    param.final$as_SD = param.SD[(1+2*no.param):(3*no.param)]
    param.final$sf_SD = param.SD[(1+3*no.param):(4*no.param)]
    # param.final$ar_SD = with(param.final, (af_SD*af_SD + as_SD*as_SD)^0.5)
    
    
    # # Calculate the parameter set from linear and quardatic equations
    # # if (no.param == 1) {
    # #   k.i = k[1]; Y.i = Y[1]; af.i = af[1]; as.i = as[1]; sf.i = sf[1]
    # # }
    # if (no.param == 2) {
    #   k.i = k[1] + k[2]*i; Y.i = Y[1]+ Y[2]*i; af.i = af[1]+ af[2]*i; as.i = as[1]+ as[2]*i; sf.i = sf[1]+ sf[2]*i
    # }
    # if (no.param == 3) {
    #   k.i = k[1] + k[2]*i + k[3]*i*i; Y.i = Y[1]+ Y[2]*i + Y[3]*i*i; af.i = af[1]+ af[2]*i + af[3]*i*i; 
    #   as.i = as[1]+ as[2]*i + as[3]*i*i; sf.i = sf[1]+ sf[2]*i + sf[3]*i*i
    # }
    
    # # Find the correlation matrix between parameter set
    # corrMatrix = cor(param.final[,c("k","Y","af","as","ar","sf")])
    # postscript(file = paste("output/figures/corrMatrix_vol_",vol[v],"_par_",no.param.par.var[z],".eps",sep=""), height=8, width=8, paper="special",
    #            family="Helvetica", fonts="Helvetica", horizontal=FALSE, onefile=FALSE)
    # corrplot(corrMatrix,tl.cex=1.5,title="Correlation Matrix", method="circle", is.corr=FALSE,type="full", cl.cex=2, 
    #          addgrid.col="blue",addshade="positive", addCoef.col = rgb(0,0,0), mar=c(0,0,1,0), diag= FALSE)
    # dev.off()
    
    
    # Calculate final output set from the predicted parameter set
    Mleaf = Mstem = Mroot = c()
    Mleaf[1] <- data$Mleaf[1]
    Mstem[1] <- data$Mstem[1]
    Mroot[1] <- data$Mroot[1]
    output.final = model(data$GPP,data$Rd,no.param,Mleaf,Mstem,Mroot,param.final$Y,
                         param.final$af,param.final$as,param.final$sf)
    
    
    # Plotting the Measured (data) vs Modelled Plant Carbon pools for plotting and comparison
    output.final$Date = data$Date
    names(output.final) = c("Mleaf.modelled","Mstem.modelled","Mroot.modelled","Date")
    melted.output = melt(output.final[,c("Mleaf.modelled","Mstem.modelled","Mroot.modelled","Date")], id.vars="Date")
    melted.data = melt(data[ , c("Mleaf","Mstem","Mroot","Date")], id.vars="Date")
    melted.data$Date = as.Date(melted.data$Date)
    melted.error = melt(data[ , c("Mleaf_SD","Mstem_SD","Mroot_SD","Date")], id.vars="Date")
    melted.error$Date = as.Date(melted.error$Date)
    melted.error$volume = as.factor(vol[v])
    melted.error$parameter = melted.data$value
    melted.error$no.param = as.factor(no.param.par.var[z])
    
    melted.output$Date = as.Date(melted.output$Date)
    melted.data$volume = as.factor(vol[v])
    melted.output$volume = as.factor(vol[v])
    melted.output$no.param = as.factor(no.param.par.var[z])
    
    
    # Calculate daily parameter values with SD
    Days <- seq(1,nrow(data), length.out=nrow(data))
    param.daily = param.final[1,] 
    
    if (no.param == 1) {
      for (i in 2:length(Days)) {
        param.daily[i,] = param.final[1,]
      }
    }
    if (no.param == 2) {
      for (i in 2:length(Days)) {
        param.daily[i,] = param.final[1,] + param.final[2,] * i
      }
    }
    # if (no.param == 3) {
    #   for (i in 2:length(Days)) {
    #     param.daily[i,] = param.final[1,] + param.final[2,] * i + param.final[3,] * i^2
    #   }
    # }
    if (no.param == 3) {
      for (i in 2:length(Days)) {
        param.daily[i,1:no.var] = param.final[1,1:no.var] + param.final[2,1:no.var] * i + param.final[3,1:no.var] * i^2
      }
    }
    for (i in (no.var+1):(2*no.var)) {
      param.daily[,i] = ((param.final[1,i]^2 + param.final[2,i]^2 + param.final[3,i]^2)/3)^0.5
    }
    param.daily$ar = 1 - param.daily$af - param.daily$as
    param.daily$ar_SD = with(param.daily, (af_SD*af_SD + as_SD*as_SD)^0.5)
    param.daily$Date = as.Date(data$Date)
    
    # Calculate the mean sf values
    # param.sf.mean$sf.mean = (param.daily$sf + param.sf.mean$sf.mean) / (v+z-1)
    
    # Find the correlation matrix between parameters and data set
    if (no.param.par.var[z] > 1) {
      # corr_data = cbind(param.daily[,c("k","Y","af","as","ar","sf")],data[,c("GPP","Rd")])
      # corrMatrix2 = cor(param.daily[,c("k","Y","af","as","ar","sf")])
      corr_data = cbind(param.daily[,c("Y","af","as","ar","sf")],data[,c("GPP","Rd","Mleaf","Mstem","Mroot")])
      corrMatrix2 = cor(corr_data, use="pairwise.complete.obs")
      
      # Create a matrix plot of scatterplots
      # pairs(param.daily[,c("k","Y","af","as","ar","sf")])
      png(height=1200, width=1200, pointsize=25, file = paste("output/figures/corr_param_data/scatterplot/param_data_",v,"_scatter_vol_",vol[v],"_par_",no.param.par.var[z],".png",sep=""))
      pairs(corr_data, use="pairwise.complete.obs",main=paste("Scatter plot of Parameter and data for vol",vol[v],"with par",no.param.par.var[z]), line.main=1.5, oma=c(2,2,3,2))
      dev.off()
      
      # Plotting the correlation matrix between parameter and data
      png(height=1200, width=1200, pointsize=25, file = paste("output/figures/corr_param_data/correlation_matrix/param_data_",v,"_corr_vol_",vol[v],"_par_",no.param.par.var[z],".png",sep=""))
      corrplot(corrMatrix2,tl.cex=1.5,title=paste("Correlation between Parameter and data for vol",vol[v],"with par",no.param.par.var[z]), method="circle", na.label = "NA", type="full", cl.cex=2,
               addgrid.col="blue",addshade="positive", addCoef.col = rgb(0,0,0), mar=c(0,0,1,0), diag= FALSE)
      dev.off()
    }
    
    # Plotting the parameter sets over time
    # param.final$Date = data$Date[seq(1,nrow(data),param.vary)]
    melted.param1 = melt(param.daily[,c("Y","af","as","ar","sf","Date")], id.vars="Date")
    melted.param2 = melt(param.daily[,c("Y_SD","af_SD","as_SD","ar_SD","sf_SD","Date")], id.vars="Date")
    melted.param = data.frame(melted.param1$Date, melted.param1$variable, melted.param1$value, melted.param2$value)
    names(melted.param) = c("Date","variable","Parameter","Parameter_SD")
    melted.param$Date = as.Date(melted.param$Date)
    melted.param$volume = vol[v]
    melted.param$no.param = as.factor(no.param.par.var[z])
    
    
    # Plotting C pools over time for individual volume and No. of parameter
    pd <- position_dodge(3) # move the overlapped errorbars horizontally
    p1 = ggplot(melted.error, aes(x=Date, y=parameter, colour=variable, group=variable)) +
      geom_errorbar(data = melted.error, aes(ymin=parameter-value, ymax=parameter+value), width=3, size=0.3) +
      geom_line(data = melted.output, aes(x = Date, y = value)) + 
      geom_point(shape = 1, size = 0.5) +
      theme_bw() +
      ylab("Plant Carbon pool (gC)") +
      ggtitle(paste("Measured (circles) vs Modelled (lines) C pools for vol",vol[v],"with par",no.param.par.var[z],"without storage")) +
      scale_colour_discrete(name="C pools",
                            breaks=c("Mleaf_SD","Mleaf.modelled", "Mroot_SD","Mroot.modelled","Mstem_SD","Mstem.modelled"),
                            labels=c("Mleaf","Mleaf.modelled", "Mroot","Mroot.modelled","Mstem","Mstem.modelled")) +
      theme(plot.title = element_text(size = 12, face = "bold")) +
      theme(legend.title = element_text(colour="chocolate", size=12, face="bold")) +
      theme(axis.title.x = element_text(size = 12, vjust=-.2)) +
      theme(axis.title.y = element_text(size = 12, vjust=0.3))
    # annotate("text", x = melted.output$Date[20], y = max(output$Mstem,na.rm = TRUE), size = 3, 
    #          label = paste("Mean k = ", round(mean(param.final[,1]), 3), "\nMean Y = ", round(mean(param.final[,2]), 3),
    #                        "\nMean af = ", round(mean(param.final[,3]), 3), "\nMean as = ", round(mean(param.final[,4]), 3),
    #                        "\nMean ar = ", round(mean(param.final[,7]), 3), "\nMean sf = ",round(mean(param.final[,5]), 3), "\nChain length = ", chainLength-bunr_in))
    p1
    ggsave(p1,filename=paste("output/figures/Cpools/Measured_vs_Modelled_Carbon_pools_",v,"_vol_",vol[v],"_par_",no.param.par.var[z],"_without_storage",".png",sep=""))
    
    
    # Plotting Allocation fractions over time for individual volume and No. of parameter
    # pd <- position_dodge(3) # move the overlapped errorbars horizontally
    p2 = ggplot(data = melted.param, aes(x = Date, y = Parameter, group = variable, colour=factor(variable))) +
      # geom_line(position=pd) +
      geom_errorbar(data = melted.param, aes(ymin=Parameter-Parameter_SD, ymax=Parameter+Parameter_SD), width=0.5, size=0.1) +
      geom_point(size=0.01) + # 21 is filled circle
      # geom_point(position=pd, size=1.5, shape=21, stroke=1.25, fill="white") + # 21 is filled circle
      xlab("Days") +
      ylab("Parameters") +
      scale_y_continuous(limits = c(-0.15,1)) +
      ggtitle(paste("Modelled allocation fractions for vol",vol[v],"with par",no.param.par.var[z],"without storage")) +
      scale_colour_hue(name="Parameter",    # Legend label, use darker colors
                       l=40) +                    # Use darker colors, lightness=40
      # scale_y_continuous(breaks=0:10*0.1)  # Set tick every 0.1
      theme_bw() +
      theme(plot.title = element_text(size = 12, face = "bold")) +
      theme(legend.title = element_text(colour="chocolate", size=12, face="bold")) +
      theme(axis.title.x = element_text(size = 12, vjust=-.2)) +
      theme(axis.title.y = element_text(size = 12, vjust=0.3))
    # + theme(legend.justification=c(1,1),
    #       legend.position=c(1,1)) # Position legend in bottom right
    p2
    ggsave(p2,filename=paste("output/figures/AF/Allocation_fractions_over_time_",v,"_vol_",vol[v],"_par_",no.param.par.var[z],"_without_storage",".png",sep=""))
    
    
    # Storing the summary of data, outputs, Cstorage, parameters
    if (q == 1) {
      summary.data = melted.data
      summary.error = melted.error
      summary.output = melted.output
      summary.param = melted.param
    }
    if (q > 1) {
      summary.output = rbind(summary.output,melted.output)
      summary.param = rbind(summary.param,melted.param)
      summary.error = rbind(summary.error,melted.error)
      if (z == 1) {
        summary.data = rbind(summary.data,melted.data)
        # summary.error = rbind(summary.error,melted.error)
      }
    }
    
    
    # Display the Acceptance rate of the chain
    nAccepted = length(unique(pChain[,1]))
    acceptance = (paste("Volume =",vol[v],", Total Parameter number =",no.param.par.var[z],": ", nAccepted, "out of ", chainLength-bunr_in, "candidates accepted ( = ",
                        round(100*nAccepted/chainLength), "%)"))
    print(acceptance)
    
    
    # # Plotting all parameter time series seperately with a moving average for the overall trend
    # ma <- function(x,n=5){filter(x,rep(1/n,n), sides=2)}
    # n = 10
    # png(file = paste("output/figures/Parameters_vol_",vol[v],"_par_",no.param.par.var[z], ".png", sep = ""))
    # par(mfrow=c(2,3))
    # plot(param.final$k,type='p',col="red",main="Utilization coefficient, k",xlab="Days")
    # # lines(ma(param.final$k,n),type='l',col="black")
    # plot(param.final$Y,type='p',col="chocolate",main="Allocation fraction to Biomass, Y",xlab="Days")
    # # lines(ma(param.final$Y,n),type='l',col="black")
    # plot(param.final$af,type='p',col="green",main="Allocation fraction to foliage, af",xlab="Days")
    # # lines(ma(param.final$af,n),type='l',col="black")
    # plot(param.final$as,type='p',col="blue",main="Allocation fraction to stem, as",xlab="Days")
    # # lines(ma(param.final$as,n),type='l',col="black")
    # plot(param.final$ar,type='p',col="magenta",main="Allocation fraction to root, ar",xlab="Days")
    # # lines(ma(param.final$ar,n),type='l',col="black")
    # plot(param.final$sf,type='p',col="purple",main="Foliage tunrover rate, sf",xlab="Days")
    # # lines(ma(param.final$sf,n),type='l',col="black")
    # dev.off()
    
    
    # Plotting all parameter whole iterations for Day 1 only to check the convergance
    png(file = paste("output/figures/Parameter_iterations/Parameter_iterations_day1_",v,"_vol_",vol[v],"_par_",no.param.par.var[z], ".png", sep = ""))
    par(mfrow=c(2,3),oma = c(0, 0, 2, 0))
    plot(pChain[,1],col="red",main="Growth respiration rate at Day 1",cex.lab = 1.5,xlab="Iterations",ylab="Y",ylim=c(param.Y[1,1],param.Y[1,3]))
    plot(pChain[,1+no.param],col="green",main="Alloc frac to foliage at Day 1",cex.lab = 1.5,xlab="Iterations",ylab="af",ylim=c(param.af[1,1],param.af[1,3]))
    plot(pChain[,1+2*no.param],col="magenta",main="Alloc frac to stem at Day 1",cex.lab = 1.5,xlab="Iterations",ylab="as",ylim=c(param.as[1,1],param.as[1,3]))
    plot(pChain[,1+3*no.param],col="blue",main="Foliage turnover at Day 1",cex.lab = 1.5,xlab="Iterations",ylab="sf",ylim=c(param.sf[1,1],param.sf[1,3]))
    plot(pChain[,1+4*no.param],col="green",main="Log-likelihood",cex.lab = 1.5,xlab="Iterations",ylab="Log-likelihood")
    title(main = paste("First day Parameter iterations for vol",vol[v],"with par",no.param.par.var[z]), outer=TRUE, cex = 1.5)
    dev.off()
    
    
    # Store the final mean parameter values
    param.mean[q,c(1:5)] = colMeans(param.daily[ , c("Y","af","as","ar","sf")])
    param.mean$volume[q] = vol[v]
    param.mean$no.param[q] = no.param.par.var[z]
    
    
    # Calcualte LogLi, AIC, BIC, Time to find the most accurate model for best balance between model fit and complexity
    output.final1 = output.final
    names(output.final1) = c("Mleaf","Mstem","Mroot","Date") # Rename for the logLikelihood function
    aic.bic[q,1] <- logLikelihood(data,output.final1) # Calculate logLikelihood
    
    k1 = 2 # k = 2 for the usual AIC
    npar = no.param*no.var # npar = total number of parameters in the fitted model
    aic.bic[q,2] = -2*aic.bic[q,1] + k1*npar
    
    n = sum(!is.na(data$Mleaf)) + sum(!is.na(data$Mstem)) + sum(!is.na(data$Mroot))
    k2 = log(n) # n being the number of observations for the so-called BIC
    aic.bic[q,3] = -2*aic.bic[q,1] + k2*npar
    
    time$end.time[q] <- Sys.time()
    time$time.taken[q] <- time$end.time[q] - time$start.time[q]
    aic.bic[q,4] = time$time.taken[q]
    aic.bic[q,5] = vol[v]
    aic.bic[q,6] = no.param.par.var[z]
  }
}

names(aic.bic) <- c("logLi","aic","bic","time","volume","no.param")
write.csv(aic.bic, file = "output/processeddata/logli_aic_bic_time_without_storage.csv", row.names = FALSE)
melted.aic.bic = melt(aic.bic, id.vars=c("no.param","volume"))

# write.csv(param.sf.mean, file = "output/processeddata/param.sf.mean.csv", row.names = FALSE)
# plot(param.daily$sf,type='l',col="red",main="Leaf turnover, sf",xlab="Days")

# This script creates the figures and saves those
source("generate_figures_CBM_v3_without_storage.R")
# 
# 
# setwd("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/output/figures/corrMatrix")
# plots1 <- lapply(ll <- list.files(patt='.*[.]png'),function(x){
#   img <- as.raster(readPNG(x))
#   rasterGrob(img, interpolate = FALSE)
# })
# ggsave("corrMatrix_multipage.pdf", marrangeGrob(grobs=plots1,nrow=2,ncol=length(no.param.par.var)))
# 
# 
# setwd("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/output/figures/AF")
# plots2 <- lapply(ll <- list.files(patt='.*[.]png'),function(x){
#   img <- as.raster(readPNG(x))
#   rasterGrob(img, interpolate = FALSE)
# })
# ggsave("AF_multipage.pdf", marrangeGrob(grobs=plots2,nrow=2,ncol=length(no.param.par.var)))
# 
# 
# setwd("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/output/figures/Cpools")
# plots3 <- lapply(ll <- list.files(patt='.*[.]png'),function(x){
#   img <- as.raster(readPNG(x))
#   rasterGrob(img, interpolate = FALSE)
# })
# ggsave("Cpools_multipage.pdf", marrangeGrob(grobs=plots3,nrow=2,ncol=length(no.param.par.var)))
# 
# 
# setwd("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/output/figures/Parameter_iterations")
# plots4 <- lapply(ll <- list.files(patt='.*[.]png'),function(x){
#   img <- as.raster(readPNG(x))
#   rasterGrob(img, interpolate = FALSE)
# })
# ggsave("Parameter_iterations_multipage.pdf", marrangeGrob(grobs=plots4,nrow=2,ncol=length(no.param.par.var)))
# 
# 
# setwd("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/output/figures/corr_param_data/correlation_matrix")
# plots5 <- lapply(ll <- list.files(patt='.*[.]png'),function(x){
#   img <- as.raster(readPNG(x))
#   rasterGrob(img, interpolate = FALSE)
# })
# ggsave("corrmatrix_param_data_multipage.pdf", marrangeGrob(grobs=plots5,nrow=2,ncol=2))
# 
# 
# setwd("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/output/figures/corr_param_data/scatterplot")
# plots6 <- lapply(ll <- list.files(patt='.*[.]png'),function(x){
#   img <- as.raster(readPNG(x))
#   rasterGrob(img, interpolate = FALSE)
# })
# ggsave("scatter_param_data_multipage.pdf", marrangeGrob(grobs=plots6,nrow=2,ncol=2))





