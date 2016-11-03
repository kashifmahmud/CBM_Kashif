# Carbon balance model 
# Developed by Kashif Mahmud and Belinda Medlyn (November 2016)
# k.mahmud@westernsydney.edu.au

# This sript reads the Pot experiment raw data to carry out Bayesian calibration for 5 variables (allocation fractions: "k","Y",af","as","sf") on 
# various temporal scales to estimate seedling Carbon pools (Cleaf,Cstem,Croot,Cstorage)
##############################

# Import daily GPP, daily Rd data
GPP.data.raw = read.csv("rawdata/GPP.csv") # Units gC d-1
Rd.data.raw = read.csv("rawdata/Rd.csv") # Units g C g-1 plant d-1
tnc.data.raw = read.csv("rawdata/tnc_fortnightly_data.csv") # Units g plant-1

# Import weekly Cleaf, weekly Cstem, initial/harvest Croot data with Mean and SD
Mleaf.data.raw = read.csv("rawdata/Cleaf_weekly_data.csv") # Units gC d-1
Mstem.data.raw = read.csv("rawdata/Cstem_weekly_data.csv") # Units gC d-1
Mroot.data.raw = read.csv("rawdata/Croot_twice_data.csv") # Units gC d-1
