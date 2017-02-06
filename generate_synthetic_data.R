# Carbon balance model 
# Developed by Kashif Mahmud and Belinda Medlyn (October 2016)
# k.mahmud@westernsydney.edu.au

# Generate synthetic data for Cleaf,Cstem,Croot,Sleaf with Mean and SD
day = 1:30 # Time in days

# Gerenate normally distributed random numbers for GPP and Rd with mean and SD
GPP.data = rnorm2(length(day),1,0.4) # Generate GPP data with mean=15, sd=3  # Units g C d-1
GPP.data = GPP.data[order(GPP.data)] # Order the genrenated data sequentially to have realistic GPP data
Rd.data = rnorm2(length(day),0.04,0.008) # Generate Rd data with mean=0.04, sd=0.008  # Units g C g-1 C d-1

# Generate random parameter set (just one value for the whole duration, no temporal variation) 
# to create synthetic Cleaf,Cstem,Croot data
k.data = 0.4 
af.data = 0.5 
as.data = 0.2
sf.data = 0.05 
Y.data = 0.3  # Allocation fraction to growth respiration (From literature: Y ~ 0.3)
param.data = c(k.data,Y.data,af.data,as.data,(1-af.data-as.data),sf.data)
  
# # Generate random parameter sets for synthetic Cleaf,Cstem,Croot data generation
# k.data = rnorm2(length(day),0.5,0.05) # Generate data 'k' values with mean=0.5, sd=0.05 
# af.data = rnorm2(length(day),1/3,0.1) # Generate data 'af' values with mean=1/3, sd=0.1 
# as.data = rnorm2(length(day),1/3,0.1) # Generate data 'as' values with mean=1/3, sd=0.1
# sf.data = rnorm2(length(day),1/50,1/500) # Generate data 'sf' values with mean=1/50, sd=1/500
# Y = 0.3  # Allocation fraction to growth respiration (Fromm literature: Y ~ 0.3)

Cstorage.data <- Mleaf.data <- Mroot.data <- Mstem.data <- c()
Cstorage.data[1] <- 0.5 # Units g C d-1
Mleaf.data[1] <- 2 # Units g C d-1
Mstem.data[1] <- 1.5 # Units g C d-1
Mroot.data[1] <- 1 # Units g C d-1

j = c()
j[1] = 0
i = seq(1,length(day),1)
# param.vary = ceiling(length(day)/no.param.par.var[z]) # How many days the parameter set remain unchanged (weekly = 7; monthly = 30; just one parameter = nrow(data))
# no.param = ceiling(length(day)/param.vary) # number of parameter set for the whole duration of experiment (121 days)
# param.vary = 121 # How many days the parameter set remain unchanged (just one parameter = nrow(data))
# no.param = ceiling(length(day)/param.vary) # number of parameter set for the whole duration of experiment (121 days)
j[i] = i - ceiling(i/121)  # j is for parameter settings for various time frames

# Generating the synthetic data sets
data = model(GPP.data,Rd.data,j,Mleaf.data,Mstem.data,Mroot.data,Y.data,k.data,af.data,as.data,sf.data)
# output.data = model(GPP.data,Rd.data,Cstorage.data,Cleaf.data,Cstem.data,Croot.data,Y,k.data,af.data,as.data,sf.data)
# names(output.data) = c("Cstorage.data","Mleaf.data","Mstem.data","Mroot.data")

# Removing some data points randomly
library(purrr)
gap.data = map_df(data, function(x) {x[sample(c(TRUE, NA), prob = c(0.2, 0.8), size = length(x), replace = TRUE)]})
gap.data[1,] = data[1,]
data = gap.data

# Initialize SD of data sets
data$Mleaf_SD = data$Mleaf * runif(30,0.05,0.1)
data$Mstem_SD = data$Mstem * runif(30,0.05,0.1)
data$Mroot_SD = data$Mroot * runif(30,0.05,0.1)
data$Sleaf_SD = data$Sleaf * runif(30,0.05,0.1)
data$Date = day
data$GPP = GPP.data
data$Rd = Rd.data

# # Plotting the synthetic data sets
# plot(data[ ,c("Mleaf","Mstem","Mroot","Sleaf")],pch=1,xlab="Days",ylab="gC",main="Different Carbon pool measurements") #plot
# legend("topleft", legend = c("Mleaf.data","Mstem.data","Mroot.data","Sleaf.data"), col=1:3, pch=0.75) # optional legend

