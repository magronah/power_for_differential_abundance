source("Load_Packages.R")
source("fitting_fun.R")
source("simulate_fun.R")
source("utils.R")
######################################################
path  =  "datasets2/"

countdata_list_obs  =  readRDS(file = paste0(path,"data.rds"))
metadata_list_obs   =  readRDS(file = paste0(path,"metadata.rds"))
######################################################
res = foreach(i = 1:7) %do% {
  source("Load_Packages.R")
  source("fitting_fun.R")
  source("utils.R")
  
  nam  =   names(metadata_list_obs)[i]
  countdata  =   countdata_list_obs[[i]]
  metadata   =   metadata_list_obs[[i]]
  ######################################################
  #filter low abundance taxa
  filter_low_otu    = filter_fun(countdata,metadata,abund_thresh=5,
                                 sample_thresh=3)
  
  countdata_filt    =  filter_low_otu$countdata_filt
  control_count     =  filter_low_otu$control_count
  treat_count       =  filter_low_otu$treat_count
  
  #run deseq for log fold change estimates 
  deseq          =   deseqfun(countdata_filt,metadata,ref_name="NT")
  countdata_norm =   deseq$normalised_count  
  dispersion     =   deseq$dispersion
  logfoldchange  =   deseq$logfoldchange
  logmean        =   deseq$logmean
  
  ## fit mean log abundance 
  mean_fit         =    logmean_fit(logmean)
  logmean_param    =    mean_fit$logmean_param
  optim_logmean    =    mean_fit$components
  logmean_density  =    mean_fit$density_plot
  
  ## fit dispersion  
  dispers_fit        =    dispersion_fit(dispersion,logmean)
  dispersion_param   =    dispers_fit$param
  dispersion_confint =    dispers_fit$confint
  
  logfoldchange_param  =    logfoldchange_fit(logmean,logfoldchange)
  
  
  list(countdata_filted =   countdata_filt,
                 countdata_norm   =   countdata_norm,
                 control_count    =   control_count,
                 treat_count      =   treat_count,
                 dispersion       =   dispersion,
                 logfoldchange    =   logfoldchange,
                 logmean          =   logmean,
                 logmean_param    =   logmean_param,
                 optim_logmean    =   optim_logmean,
                 logmean_density  =   logmean_density,
                 dispersion_param =   dispersion_param,
                 dispersion_confint   = dispersion_confint,
                 logfoldchange_param  =   logfoldchange_param)
}
##################################################################
names(res)        =   names(countdata_list_obs)

otu_dataset_list  = list(filtered_otu    =  read_data(res, "countdata_filted"),
                         normalised_otu  =  read_data(res, "countdata_norm"),
                         control_otu     =  read_data(res, "control_count"),
                         treatment_otu   =  read_data(res, "treat_count"))

estimates_obs_data_list =  list(logfoldchange  =  read_data(res, "logfoldchange"),
                                dispersion     =  read_data(res, "dispersion"),
                                logmean        =  read_data(res, "logmean"))

fit_param_list  =  list(logmean_param      =  read_data(res, "logmean_param"),
                        optim_logmean        =  read_data(res, "optim_logmean"),
                        dispersion_param     =  read_data(res, "dispersion_param"),
                        dispersion_confint   =  read_data(res, "dispersion_confint"),
                        logfoldchange_param  =  read_data(res, "logfoldchange_param"))

saveRDS(otu_dataset_list, file = paste0(path,"otu_dataset_list.rds"))
saveRDS(estimates_obs_data_list, file = paste0(path,"estimates_obs_data_list.rds"))
saveRDS(fit_param_list, file = paste0(path,"fit_param_list.rds"))

