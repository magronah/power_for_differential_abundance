######################################################
##' This script reads in all the results from save_estimates.R
######################################################
path = paste0(getwd(),"/datasets/")
###########################Read otu dataset###########################
   otu_dataset_list     =    readRDS(paste0(path,"otu_dataset_list.rds"))
   filtered_otu_list    =    otu_dataset_list[["filtered_otu"]]
  normalised_otu_list   =    otu_dataset_list[["normalised_otu"]]
     control_otu_list   =    otu_dataset_list[["control_otu"]]
   treatment_otu_list   =    otu_dataset_list[["treatment_otu"]]

##########Read logfoldchange, logmean and dispersion estimates#########
estimates_obs_data_list  =   readRDS( paste0(path,"estimates_obs_data_list.rds"))
logfoldchange_list       =   estimates_obs_data_list[["logfoldchange"]]
logmean_list             =   estimates_obs_data_list[["logmean"]]
dispersion_list          =   estimates_obs_data_list[["dispersion"]]

########################Read parameters################################
fit_param_list           =   readRDS(paste0(path,"fit_param_list.rds"))
logmean_param_list       =   fit_param_list[["logmean_param"]]
dispersion_param_list    =   fit_param_list[["dispersion_param"]]
logfoldchange_param_list =   fit_param_list[["logfoldchange_param"]]







