# Carbon balance model (CBM)
# Developed by Kashif Mahmud and Belinda Medlyn (November 2016)
# k.mahmud@westernsydney.edu.au

# This code carries out Bayesian calibration for 4 variables (allocation fractions: "k","af","as","sf") on 
# daily time scale (e.g. 120 days) to estimate Carbon pools (Cstorage,Cleaf,Cstem,Croot)

##############################
# Version = v15: MCMC with soil manipulation pot experiment data for all treatments (including the free seedling), 
# This version considers either daily/weekly/monthly/just one parameter set for 5 variables ("k","Y","af","as","sf")
# So we can set the parameters for various time frames
# Also calculates the MCMC SDs for all parameters with different time frames, and also the LogLi, AIC, BIC, time (measures 
# for best model selection) to select the best parameter set
# Finally save the figures in Github/results folder
##############################
# Set working directory for saving figures
setwd("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif")

# This script cleans the workspace, loads necessary Rfunctions and packages
source("MCMC_CBM_load.R")

# This sript reads the Pot experiment raw data
source("MCMC_CBM_readdata.R")

# Assign inputs for MCMC
chainLength = 1500 # Setting the length of the Markov Chain to be generated
no.var = 5 # variables to be modelled are: k,Y,af,as,sf

# Assign pot volumes and number of parameters per varible in temporal scale
vol = c(1000) # test run
no.param.par.var = c(3) # test run
# GPP.data.raw = read.csv("GPP.csv") # Units gC d-1
# vol = unique(GPP.data.raw$volume)[order(unique(GPP.data.raw$volume))]
# no.param.par.var = c(1,2,3,4,5,6,9) # temporal parameter count per variable

param.mean = data.frame(matrix(ncol = no.var+1, nrow = length(no.param.par.var)*length(vol)))
names(param.mean) = c("k","Y","af","as","ar","sf")
aic.bic = data.frame(matrix(ncol = 4, nrow = length(no.param.par.var)*length(vol)))
time = data.frame(no.param=rep(no.param.par.var,length(vol)),
                  start.time=numeric(length(no.param.par.var)*length(vol)),
                  end.time=numeric(length(no.param.par.var)*length(vol)),
                  time.taken=numeric(length(no.param.par.var)*length(vol)))
q = 0 # Indicates the iteration number


for (z in 1:length(no.param.par.var)) {
  for (v in 1:length(vol)) {
    # This script process the raw data
    source("MCMC_CBM_processdata.R")
    
    
    # Initialize few output data files
    q = q + 1
    time$start.time[q] <- Sys.time()
    param.vary = ceiling(nrow(data)/no.param.par.var[z]) # How many days the parameter set remain unchanged (weekly = 7; monthly = 30; just one parameter = nrow(data))
    no.param = ceiling(nrow(data)/param.vary) # number of parameter set for the whole duration of experiment (121 days)
    j = c()
    j[1] = 0
    i = seq(1,nrow(data),1)
    j[i] = i - ceiling(i/param.vary)*1  # j is for parameter settings for various time frames
    
    
    # Setting lower and upper bounds of the prior parameter pdf, and starting point of the chain
    param.k <- matrix(c(0,0.45,1) , nrow=no.param, ncol=3, byrow=T) 
    param.Y <- matrix(c(0.1,0.3,0.5) , nrow=no.param, ncol=3, byrow=T) 
    param.af <- matrix(c(0,0.45,0.7) , nrow=no.param, ncol=3, byrow=T) 
    param.as <- matrix(c(0,0.2,0.5) , nrow=no.param, ncol=3, byrow=T) 
    param.sf <- matrix(c(0,0.05,0.1) , nrow=no.param, ncol=3, byrow=T) 
    
    param = data.frame(param.k,param.Y,param.af,param.as,param.sf)
    names(param) <- c("k_min","k","k_max","Y_min","Y","Y_max","af_min","af","af_max","as_min","as","as_max","sf_min","sf","sf_max")
    pMinima <- param[ ,c("k_min","Y_min","af_min","as_min","sf_min")]
    pMaxima <- param[ ,c("k_max","Y_max","af_max","as_max","sf_max")]
    pValues <- param[ ,c("k","Y","af","as","sf")] # Starting point of the chain
    pChain <- matrix(0, nrow=chainLength, ncol = no.param*no.var+1) # Initialising the chain
    
    
    # Defining the variance-covariance matrix for proposal generation
    vcov = (0.01*(pMaxima-pMinima))^2
    vcovProposal =  vcov[1,] # The higher the coefficient, the lower the acceptance rate with better matching
    
    
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
    output = model(data$GPP,data$Rd,j,Mleaf,Mstem,Mroot,pValues$Y,pValues$k,pValues$af,pValues$as,pValues$sf)
    logL0 <- logLikelihood(data,output) # Calculate log likelihood of starting point of the chain
    pChain[1,] <- c(pValues$k,pValues$Y,pValues$af,pValues$as,pValues$sf,logL0) # Assign the first parameter set with log likelihood
    
    
    # Calculating the next candidate parameter vector, as a multivariate normal jump away from the current point
    for (c in (2 : chainLength)) {
      candidatepValues = matrix(ncol = no.var, nrow = no.param)
      for (i in 1:no.var) {
        candidatepValues[ , i] = rmvnorm(n=1, mean=pValues[ , i],
                                         sigma=diag(vcovProposal[i],no.param)) 
      }
      candidatepValues = data.frame(candidatepValues)
      names(candidatepValues) <- c("k","Y","af","as","sf")
      
      
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
        out.cand = model(data$GPP,data$Rd,j,Mleaf,Mstem,Mroot,candidatepValues$Y,
                         candidatepValues$k,candidatepValues$af,candidatepValues$as,candidatepValues$sf)
        logL1 <- logLikelihood(data,out.cand) # Calculate log likelihood

        
        # Calculating the logarithm of the Metropolis ratio
        logalpha <- (logPrior1+logL1) - (logPrior0+logL0) 
        # Accepting or rejecting the candidate vector
        if ( log(runif(1, min = 0, max =1)) < logalpha ) { 
          pValues <- candidatepValues
          logPrior0 <- logPrior1
          logL0 <- logL1
        }
      }
      pChain[c,] <- c(pValues$k,pValues$Y,pValues$af,pValues$as,pValues$sf,logL0)
    }
    # Discard the first 500 iterations for Burn-IN in MCMC
    pChain <- pChain[501:nrow(pChain),]
    
    # Store the final parameter set values
    param.set = colMeans(pChain[ , 1:(no.param*no.var)])
    param.SD = apply(pChain[ , 1:(no.param*no.var)], 2, sd)
    param.final = data.frame(matrix(ncol = (no.var+1)*2, nrow = no.param))
    names(param.final) <- c("k","Y","af","as","ar","sf","k_SD","Y_SD","af_SD","as_SD","ar_SD","sf_SD")
    param.final$k = param.set[1:no.param]
    param.final$Y = param.set[(1+no.param):(2*no.param)]
    param.final$af = param.set[(1+2*no.param):(3*no.param)]
    param.final$as = param.set[(1+3*no.param):(4*no.param)]
    param.final$sf = param.set[(1+4*no.param):(5*no.param)]
    param.final$ar = 1 - param.final$af - param.final$as
    param.final$k_SD = param.SD[1:no.param]
    param.final$Y_SD = param.SD[(1+no.param):(2*no.param)]
    param.final$af_SD = param.SD[(1+2*no.param):(3*no.param)]
    param.final$as_SD = param.SD[(1+3*no.param):(4*no.param)]
    param.final$sf_SD = param.SD[(1+4*no.param):(5*no.param)]
    param.final$ar_SD = with(param.final, (af_SD*af_SD + as_SD*as_SD)^0.5)
    
    # Calculate final output set from the predicted parameter set
    Mleaf = Mstem = Mroot = c()
    Mleaf[1] <- data$Mleaf[1]
    Mstem[1] <- data$Mstem[1]
    Mroot[1] <- data$Mroot[1]
    output.final = model(data$GPP,data$Rd,j,Mleaf,Mstem,Mroot,param.final$Y,
                         param.final$k,param.final$af,param.final$as,param.final$sf)
    
    
    # Plotting the Measured (data) vs Modelled Plant Carbon pools for plotting and comparison
    output.final$Date = data$Date
    names(output.final) = c("Cstorage.modelled","Mleaf.modelled","Mstem.modelled","Mroot.modelled","Sleaf.modelled","Date")
    melted.output = melt(output.final[,c("Mleaf.modelled","Mstem.modelled","Mroot.modelled","Sleaf.modelled","Date")], id.vars="Date")
    melted.Cstorage = output.final[,c("Cstorage.modelled","Date")]
    melted.data = melt(data[ , c("Mleaf","Mstem","Mroot","Sleaf","Date")], id.vars="Date")
    melted.data$Date = as.Date(melted.data$Date)
    melted.error = melt(data[ , c("Mleaf_SD","Mstem_SD","Mroot_SD","Sleaf_SD","Date")], id.vars="Date")
    melted.error$Date = as.Date(melted.error$Date)
    melted.error$volume = as.factor(vol[v])
    melted.error$parameter = melted.data$value
    melted.output$Date = as.Date(melted.output$Date)
    melted.data$volume = as.factor(vol[v])
    melted.output$volume = as.factor(vol[v])
    melted.output$no.param = as.factor(no.param.par.var[z])
    melted.Cstorage$Date = as.Date(melted.Cstorage$Date)
    melted.Cstorage$volume = as.factor(vol[v])
    melted.Cstorage$no.param = as.factor(no.param.par.var[z])
    
    
    # Plotting the parameter sets over time
    param.final$Date = data$Date[seq(1,nrow(data),param.vary)]
    melted.param1 = melt(param.final[,c("k","Y","af","as","ar","sf","Date")], id.vars="Date")
    melted.param2 = melt(param.final[,c("k_SD","Y_SD","af_SD","as_SD","ar_SD","sf_SD","Date")], id.vars="Date")
    melted.param = data.frame(melted.param1$Date, melted.param1$variable, melted.param1$value, melted.param2$value)
    names(melted.param) = c("Date","variable","Parameter","Parameter_SD")
    melted.param$Date = as.Date(melted.param$Date)
    melted.param$volume = vol[v]
    melted.param$no.param = as.factor(no.param.par.var[z])

    
    # Plotting C pools over time for individual volume and No. of parameter
    pd <- position_dodge(3) # move the overlapped errorbars horizontally
    p1 = ggplot(melted.error, aes(x=Date, y=parameter, colour=variable, group=variable)) + 
      geom_errorbar(aes(ymin=parameter-value, ymax=parameter+value,colour=variable, linetype=variable), width=7) +
      geom_line(data = melted.output, aes(x = Date, y = value, group = variable, colour=variable)) + 
      geom_point(shape = 1, size = 1, stroke = 1.25) +
      ylab("Plant Carbon pool (gC)") +
      ggtitle("Measured (circles) vs Modelled (lines) Plant Carbon pools") +
      labs(linetype="Data uncertainty",colour="C pools") +
      theme(legend.title = element_text(colour="chocolate", size=10, face="bold")) +
      annotate("text", x = melted.output$Date[20], y = max(output$Mstem,na.rm = TRUE), size = 3, 
                          label = paste("Mean k = ", round(mean(param.final[,1]), 3), "\nMean Y = ", round(mean(param.final[,2]), 3),
                                        "\nMean af = ", round(mean(param.final[,3]), 3), "\nMean as = ", round(mean(param.final[,4]), 3),
                                        "\nMean ar = ", round(mean(param.final[,7]), 3), "\nMean sf = ",round(mean(param.final[,5]), 3), "\nChain length = ", chainLength-500))
    p1
    ggsave(p1,filename=paste("output/figures/Measured_vs_Modelled_Carbon_pools_vol_",vol[v],"_par_",no.param.par.var[z],".png",sep=""))
    
    
    # Plotting Allocation fractions over time for individual volume and No. of parameter
    pd <- position_dodge(3) # move the overlapped errorbars horizontally
    p2 = ggplot(data = melted.param, aes(x = Date, y = Parameter, group = variable, colour=factor(variable))) +
      geom_line(position=pd) +
      geom_errorbar(data = melted.param, aes(ymin=Parameter-Parameter_SD, ymax=Parameter+Parameter_SD), width=5, position=pd) +
      geom_point(position=pd, size=1.5, shape=21, stroke=1.25, fill="white") + # 21 is filled circle
      xlab("Days") +
      ylab("Parameters") +
      ggtitle("Modelled allocation fractions") +
      scale_colour_hue(name="Parameter",    # Legend label, use darker colors
                       l=40) +                    # Use darker colors, lightness=40
      scale_y_continuous(breaks=0:10*0.1)  # Set tick every 0.1
    theme_bw() +
      theme(legend.justification=c(1,1),
            legend.position=c(1.1,1.1)) # Position legend in bottom right
    p2
    ggsave(p2,filename=paste("output/figures/Allocation_fractions_over_time_vol_",vol[v],"_par_",no.param.par.var[z],".png",sep=""))
    
    # Storing the summary of data, outputs, Cstorage, parameters
    if (q == 1) {
      summary.data = melted.data
      summary.error = melted.error
      summary.output = melted.output
      summary.Cstorage = melted.Cstorage
      summary.param = melted.param
    }
    if (q > 1) {
      summary.output = rbind(summary.output,melted.output)
      summary.Cstorage = rbind(summary.Cstorage,melted.Cstorage)
      summary.param = rbind(summary.param,melted.param)
      if (z == 1) {
        summary.data = rbind(summary.data,melted.data)
        summary.error = rbind(summary.error,melted.error)
      }
    }
    
    
    # Display the Acceptance rate of the chain
    nAccepted = length(unique(pChain[,1]))
    acceptance = (paste("Volume =",vol[v],", Total Parameter number =",no.param.par.var[z],": ", nAccepted, "out of ", chainLength, "candidates accepted ( = ",
                        round(100*nAccepted/chainLength), "%)"))
    print(acceptance)
    
    
    # Plotting all parameter time series seperately with a moving average for the overall trend
    ma <- function(x,n=5){filter(x,rep(1/n,n), sides=2)}
    n = 10
    png(file = paste("output/figures/Parameters_vol_",vol[v],"_par_",no.param.par.var[z], ".png", sep = ""))
    par(mfrow=c(2,3))
    plot(param.final$k,type='p',col="red",main="Utilization coefficient, k",xlab="Days")
    # lines(ma(param.final$k,n),type='l',col="black")
    plot(param.final$Y,type='p',col="chocolate",main="Allocation fraction to Biomass, Y",xlab="Days")
    # lines(ma(param.final$Y,n),type='l',col="black")
    plot(param.final$af,type='p',col="green",main="Allocation fraction to foliage, af",xlab="Days")
    # lines(ma(param.final$af,n),type='l',col="black")
    plot(param.final$as,type='p',col="blue",main="Allocation fraction to stem, as",xlab="Days")
    # lines(ma(param.final$as,n),type='l',col="black")
    plot(param.final$ar,type='p',col="magenta",main="Allocation fraction to root, ar",xlab="Days")
    # lines(ma(param.final$ar,n),type='l',col="black")
    plot(param.final$sf,type='p',col="purple",main="Foliage tunrover rate, sf",xlab="Days")
    # lines(ma(param.final$sf,n),type='l',col="black")
    dev.off()
    
    
    # Plotting all parameter whole iterations for Day 1 only to check the convergance
    png(file = paste("output/figures/Parameter_iterations_day1_vol_",vol[v],"_par_",no.param.par.var[z], ".png", sep = ""))
    par(mfrow=c(2,3))
    plot(pChain[,1],col="red",main="Utilization coefficient at Day 1",xlab="Iterations",ylab="k")
    plot(pChain[,1+no.param],col="green",main="Alloc frac to Biomass at Day 1",xlab="Iterations",ylab="Y")
    plot(pChain[,1+2*no.param],col="magenta",main="Alloc frac to foliage at Day 1",xlab="Iterations",ylab="af")
    plot(pChain[,1+3*no.param],col="blue",main="Alloc frac to stem at Day 1",xlab="Iterations",ylab="as")
    plot(pChain[,1+4*no.param],col="green",main="Foliage turnover at Day 1",xlab="Iterations",ylab="sf")
    plot(pChain[,1+5*no.param],col="magenta",main="Log-likelihood",xlab="Iterations",ylab="Log-likelihood")
    dev.off()
    
    
    # Store the final mean parameter values
    param.mean[q,c(1:6)] = colMeans(param.final[ , c(1:6)])
    param.mean$volume[q] = vol[v]
    param.mean$no.param[q] = no.param.par.var[z]
    
    
    # Calcualte LogLi, AIC, BIC, Time to find the most accurate model for best balance between model fit and complexity
    output.final1 = output.final
    names(output.final1) = c("Cstorage","Mleaf","Mstem","Mroot","Sleaf","Date") # Rename for the logLikelihood function
    aic.bic[q,1] <- logLikelihood(data,output.final1) # Calculate logLikelihood
    
    k1 = 2 # k = 2 for the usual AIC
    npar = no.param*no.var # npar = total number of parameters in the fitted model
    aic.bic[q,2] = -2*aic.bic[q,1] + k1*npar
    
    n = sum(!is.na(data$Sleaf)) + sum(!is.na(data$Mleaf)) + sum(!is.na(data$Mstem)) + sum(!is.na(data$Mroot))
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
write.csv(aic.bic, file = "output/processeddata/logli_aic_bic_time.csv", row.names = FALSE)
melted.aic.bic = melt(aic.bic, id.vars=c("no.param","volume"))


# This script creates the figures and saves those
source("MCMC_CBM_figures.R")


