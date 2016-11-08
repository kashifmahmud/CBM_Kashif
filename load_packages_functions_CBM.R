# Carbon balance model 
# Developed by Kashif Mahmud and Belinda Medlyn (November 2016)
# k.mahmud@westernsydney.edu.au

# This script cleans the workspace, loads necessary Rfunctions and packages
##############################
rm(list=ls()) # Cleaning the workspace

# Load the function to define the CBM equations to iteratively calculate Cstorage, Cleaf, Cstem, Croot, Sleaf, Sstem, Sroot
source("Rfunctions/CBM_model.R")

# Load the function to calcualte LogLi to find the most accurate model
source("Rfunctions/logLikelihood.R")

# Load the function to calcualte LogLi to find the most accurate model
source("Rfunctions/multiplot.R")

# Load the function to calcualte LogLi to find the most accurate model
source("Rfunctions/rnorm2.R")

# Load packages
# install.packages("mvtnorm")
library("mvtnorm") # Creates candidate parameter vector as a multivariate normal jump away from the current candidate
library("reshape2")
library("ggplot2")
library("lubridate")
library("rio")
library("dplyr")
library("zoo")
library("doBy")
library("corrplot")
library("png")
library("grid")
library(gridExtra)
