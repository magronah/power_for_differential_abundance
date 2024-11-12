library(foreach)
library(doParallel)
library(ggokabeito)
path = "datasets2/"
####################################################################
countdata_list_obs  =  readRDS(file = paste0(path,"data.rds"))
metadata_list_obs   =  readRDS(file = paste0(path,"metadata.rds"))
######################################################
ncores  =  10
cl      =   makeCluster(ncores)
registerDoParallel(cl)
hmp = metasparsim = list()
nsim = 1;     disp_scale = 0.3 # scaled dispersion value

res <-  foreach(n = 1:7, .errorhandling = "pass", 
                .packages = c("HMP","DirichletMultinomial","metaSPARSim")) %dopar% {  

  source("read_estimate_results.R")
  source("Load_Packages.R")
  source("simulate_fun.R")
  source("utils.R")
  
  nam    =   names(metadata_list_obs)[n]
  countdata   =   filtered_otu_list[[n]]
  metadata    =   metadata_list_obs[[n]]
  
  logmean_param = logmean_param_list[[n]]
  logfoldchange_param = logfoldchange_param_list[[n]]
  dispersion_param  =  dispersion_param_list[[n]]
  ######################################################
  #HMP needs the data sets to have sample by taxa format.
  countdata_hmp =    t(countdata)
  shape         =    dirmult(countdata_hmp)$gamma
  size          =    rowSums(countdata_hmp)
  sim_hmp       =    Dirichlet.multinomial(size, shape)
  sim_hmp       =    as.data.frame(t(sim_hmp))
  hmp[[n]]      =    sim_hmp
  # ######################################################
  # #simulate from metaSPARSim
  # #metaSPARSim needs the data sets to have taxa by sample format.
  normalised_data =    normalised_otu_list[[n]]
  control_otu     =    control_otu_list[[n]]
  treatment_otu   =    treatment_otu_list[[n]]
  nt_index        =    which(colnames(normalised_data) %in%
                               colnames(control_otu))
  asd_index       =    which(colnames(normalised_data) %in%
                               colnames(treatment_otu))
  stopifnot(colnames(control_otu) ==
              colnames(normalised_data)[nt_index])
  stopifnot(colnames(treatment_otu) ==
              colnames(normalised_data)[asd_index])

  grp_index      =    list(NT = nt_index, ASD = asd_index)
  params         =    estimate_parameter_from_data(countdata, normalised_data, grp_index)
  names(params)  =    names(grp_index)
  sim_data       =    metaSPARSim(params)

  countdata_metaSPARSim = sim_data$counts
  metasparsim[[n]] =   countdata_metaSPARSim
  #################################################################
  #########my simulating from our method 
  notu = nrow(countdata)
  ncont = length(nt_index); ntreat = length(asd_index)

  otu_data   =  countdata_sim_fun(logmean_param, logfoldchange_param, 
                                   dispersion_param,
                                   nsamp_per_group=NULL,
                                   disp_scale = disp_scale,
                                   ncont = ncont,ntreat = ntreat,
                                   notu, nsim, 
                                   seed = 100)
  
  countdata_scaled =  otu_data$countdata_list[[1]]
  
  ##########################################################
  countdata_sim_list  =  list(metaSPARSim   =   countdata_metaSPARSim,
                              our_method    =   countdata_scaled,
                              HMP           =   sim_hmp)
  
  countdata_sim_list
}

stopCluster(cl)
names(res)  =  names(metadata_list_obs) 
saveRDS(res, file = paste0(path,"HMP_metaSPARSim_sim_compare.rds"))

HMP_list           =    read_data(res, "HMP")
metaSPARSim_list   =    read_data(res, "metaSPARSim")
our_method_list    =    read_data(res, "our_method")

### Generate comparison plot
plt =  foreach(i = 1:length(res), .errorhandling = "pass") %do% {
                  nam       =    names(filtered_otu_list)[i]
          countdata_obs     =    filtered_otu_list[[i]]
                  HMP       =    HMP_list[[i]]
           metaSPARSim      =    metaSPARSim_list[[i]]
           our_method       =    our_method_list[[i]] 

                      
      
        countdata_sim       =    list(HMP = HMP, metaSPARSim = metaSPARSim,
                                     our_method =  our_method) 
                                     
  okabe_ito_colors = c("#556B2F", "#E23D28", "#0000FF","#000000")
  p1   =   compare_dataset(countdata_sim,countdata_obs,method = "var") 
  p1   =   p1  + ggtitle(nam) + scale_color_manual(values = okabe_ito_colors)  

  p2   =   compare_dataset(countdata_sim,countdata_obs,method = "mean")
  p2   =   p2  + ggtitle(nam) + scale_color_manual(values = okabe_ito_colors)  

  pp   =   list(var_plt  =  p1,  mean_plt = p2)
  
  pp
}

## Extract plots for mean and variance distributions
mean_plt  =  read_data(plt, "mean_plt")
var_plt   =  read_data(plt, "var_plt")

## show plots
pp1  =   (mean_plt[[1]]|mean_plt[[2]]|mean_plt[[3]]) +  plot_layout(guides = "collect")  
pp2  =   (var_plt[[1]]|var_plt[[2]]|var_plt[[3]])   +  plot_layout(guides = "collect")  

pp1/pp2

 


