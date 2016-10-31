CBM_Kashif

AUTHOR: Kashif Mahmud and Belinda Medlyn
k.mahmud@westernsydney.edu.au
01/11/2016

Carbon Balance Model (CBM)

A simple Carbon Balance Model (CBM) for an individual free growing plant seedling (Eucalyptus Tereticornis) for daily time scale:

The repository includes: 
i) a R project "CBM_Kashif.Rproj" to perform MCMC with a carbon balance model (CBM). 
ii) a R script "MCMC_CBM_analysis.R" to perform MCMC with a carbon balance model (CBM). 
iii) other R scripts to perform different portions of the analysis and are linked with the main script "MCMC_CBM_analysis.R".

"Instruction_files" provides a brief instructions to describe the Carbon Balance Model (CBM) structure, associated equations and the steps the R scripts are performing to process the necessary data. Another instruction file describes the MCMC parameter settings and the synthetic data generation.

All necessary raw data files that are called in the R scripts are stored in "rawdata" folder for easy access.

The "outcome" folder holds all the results (including figures, processed data tables) produced from the CBM MCMC anlysis. 

All the functions used in the R project are stored in the "Rfunctions" folder.
