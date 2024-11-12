source("reproducible/power/Load_Packages.R")
source("reproducible/power/utils.R")
source("reproducible/power/read_dataset.R")
source("reproducible/power/simulate_fun.R")
path  =  "reproducible/power/datasets2/"
############################################################
metadata_list_obs   =  readRDS(file = paste0(path,"metadata.rds"))

nsim = 2;   notu   =  1000; nsamp  =  100
disp_vec = c(0.1,0.3,0.9)
############################################################
dd_sig  =  dd_data = list()

for(i in 1:7){
                nam       =   names(metadata_list_obs)[i]
  dispersion_param      =   dispersion_param_list[[i]]
  logmean_param         =   logmean_param_list[[i]]
  logfoldchange_param   =   logfoldchange_param_list[[i]]
  countdata_filt        =   filtered_otu_list[[i]]
  
  for(j in 1:length(disp_vec)){
    disp_scale = disp_vec[j]
    
    ###### Data simulation 
    otu_data   =  countdata_sim_fun(logmean_param, logfoldchange_param, 
                                    dispersion_param,
                                    nsamp_per_group=nsamp,
                                    disp_scale = disp_scale,
                                    ncont = NULL,ntreat = NULL,
                                    notu, nsim, 
                                    seed = 100)
    
    
    countdata_list        =   otu_data$countdata_list
    metadata_list         =   otu_data$metadata_list
  
    pp1 = lapply(countdata_list, function(x){apply(x,1,var)})
    pp2 = lapply(countdata_list, function(x){apply(x,1,mean)})
    
   dd_sig[[j]]    =  data.frame(var  =  log(unlist(pp1)), 
                        mean  =  log(unlist(pp2)),
                       scale  =  rep(disp_scale, length(unlist(pp1))))
  }
   DD    =   do.call(rbind,dd_sig)
   DD$names  =  rep(nam, nrow(DD))
   dd_data[[i]] =  DD
}


df = do.call(rbind,dd_data)
df$scale  = as.factor(df$scale)


dd_list <- Map(function(x, name) {
  dd <- data.frame(mean = log(apply(x, 1, mean)), var = log(apply(x, 1, var)))
  dd$names <- rep(name, nrow(dd))
  dd
}, filtered_otu_list, names(filtered_otu_list))

dd_lis =  do.call(rbind, dd_list)

ddf = df[(df$scale) == "0.3",]

n  = 10
p1 <- ggplot(ddf, aes(x = var, colour = scale)) + 
  geom_density(lwd = 1) +
  geom_density(data = dd_lis, aes(x = var), 
               colour = "black", linetype = "dashed", lwd = 2) +
  theme_bw(base_size = n) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = n, family = "Roboto"),
    axis.text.x = element_text(family = "Roboto", size = n, color = "black"),
    axis.text.y = element_text(family = "Roboto", size = n, color = "black")
  ) +
  facet_wrap(~names, scales = "free")
p1

p2 <- ggplot(ddf, aes(x = mean, colour = scale)) + 
  geom_density(lwd = 0.8) +
  geom_density(data = dd_lis, aes(x = mean), 
               colour = "black", linetype = "dashed", lwd = 2) +
  theme_bw(base_size = n) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = n, family = "Roboto"),
    axis.text.x = element_text(family = "Roboto", size = n, color = "black"),
    axis.text.y = element_text(family = "Roboto", size = n, color = "black")
  ) +
  facet_wrap(~names, scales = "free")

p2
####################################################

### compare variance of simulation and observed
plt1[[i]]  = compare_dataset(countdata_sim_list =  countdata_list[1:5],
                             countdata_filt, method = "var",n=11)

plt2[[i]] = compare_dataset(countdata_sim_list =  countdata_list[1:5],
                            countdata_filt, method = "mean",n=11)



plt1 =  plt2 =  dd_sig = list()
disp_scale = 0.3; nsim = 3
notu   =  1000; nsamp  =  100
for(i in 1:7){
  nam  =   names(metadata_list_obs)[i]
  
  dispersion_param      =   dispersion_param_list[[i]]
  logmean_param         =   logmean_param_list[[i]]
  logfoldchange_param   =   logfoldchange_param_list[[i]]
  countdata_filt        =   filtered_otu_list[[i]]
  
    ###### Data simulation 
    otu_data   =  countdata_sim_fun(logmean_param, logfoldchange_param, 
                                    dispersion_param,
                                    nsamp_per_group=nsamp,
                                    disp_scale = disp_scale,
                                    ncont = NULL,ntreat = NULL,
                                    notu, nsim, 
                                    seed = 100)
    
    
    countdata_list      =   otu_data$countdata_list
    metadata_list         =   otu_data$metadata_list
    true_lmean_list       =   otu_data$logmean_list
    true_lfoldchange_list    =   otu_data$logfoldchange_list
    treat_countdata_list    =   otu_data$treat_countdata_list
    control_countdata_list   =   otu_data$control_countdata_list
    
    pp1 = lapply(countdata_list, function(x){apply(x,1,var)})
    pp2 = lapply(countdata_list, function(x){apply(x,1,mean)})
   
    ### compare variance of simulation and observed
    plt1[[i]]  = compare_dataset(countdata_sim_list =  countdata_list[1:nsim],
                                 countdata_filt, method = "var",n=11)
    
    plt2[[i]] = compare_dataset(countdata_sim_list =  countdata_list[1:nsim],
                                countdata_filt, method = "mean",n=11)
}

(plt1[[1]] | plt1[[2]] | plt1[[3]]) /
  (plt1[[4]] | plt1[[5]] | plt1[[6]]) /
  (plt1[[7]] | plot_spacer() | plot_spacer()) + 
  plot_layout(guides = "collect") 


(plt2[[1]] | plt2[[2]] | plt2[[3]]) /
  (plt2[[4]] | plt2[[5]] | plt2[[6]]) /
  (plt2[[7]] | plot_spacer() | plot_spacer()) + 
  plot_layout(guides = "collect") 


CV = disp = list()
res = foreach(i = 1:7, .combine = "rbind") %do%{
  nam         =   names(logmean_list)[i]
  mean        =   2^logmean_list[[i]]
  dispersion  =   dispersion_list[[i]]
  CV[[i]]     =   data.frame(cv =  coeff_var(mean, dispersion), 
                             data_accession_no  =   rep(nam, length(dispersion)))
  disp[[i]]   =   data.frame(dispersion = dispersion, 
                             data_accession_no  =   rep(nam, length(dispersion)))
}


CV = do.call(rbind,CV)
disp = do.call(rbind,disp)
ggplot(CV, aes(x = data_accession_no, y= cv)) +
  geom_violin(scale = "width") +
  geom_point(alpha= 0.1) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size =10, family = "Roboto"),
    axis.text.x = element_text(family = "Roboto",angle = 90, size =10, color = "black"),
    axis.text.y = element_text(family = "Roboto", size =10, color = "black")
  )+ 
  labs(x="data accession number", y ="coefficient of variation")



ggplot(CV, aes(x = data_accession_no, y= cv)) +
  geom_violin(scale = "width") +
  geom_point(alpha= 0.1) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size =10, family = "Roboto"),
    axis.text.x = element_text(family = "Roboto",angle = 90, size =10, color = "black"),
    axis.text.y = element_text(family = "Roboto", size =10, color = "black")
  )+ 
  labs(x="data accession number", y ="coefficient of variation")


####################################################################################
#Fit lm
combined_list <- mapply(
  function(lmc, lfc) list(lmc = lmc, lfc = lfc),
  lmc = logmean_list,
  lfc = logfoldchange_list,
  SIMPLIFY = FALSE
)

mod=  lapply(combined_list, function(x){lm(x$lfc ~ x$lmc)})

summ_mod = lapply(mod, function(x){summary(x)})
names(summ_mod)

summ_mod[["PRJNA355023"]]
summ_mod[["PRJNA453621"]]
summ_mod[["PRJNA644763"]]

##########################################################################
###############Deseq's estimates for simulated data####################
#nsamp_vec     =    seq(50,300,50)

# nsamp_vec     =   seq(10,100,10)
# contour_data  =   power_estimate =  list()
# df1   =   df2 =   list()
# 
# folder_path  =   paste0(path,"contour_data/")
# 
# for(i in 1:7){
#   nam     =   names(filtered_otu_list)[i]
# 
#   for(j in 1:length(nsamp_vec)){
#    nsam    =   nsamp_vec[j]
#    dd1      =   readRDS(paste0(folder_path,nsam,"_samples/combined_data_",
#                        nam, "_.rds"))
#    dd2      =   readRDS(paste0(folder_path,nsam,"_samples/power_estimate_",
#                                nam, "_.rds"))
# 
#    dd1$sample_size = rep(nsam, nrow(dd1))
#    dd2$sample_size = rep(nsam, nrow(dd2))
# 
#    df1[[j]]   =  dd1
#    df2[[j]]   =  dd2
#   }
# 
#   names(df1)  =  names(df2)  =  paste0(nsamp_vec, "samples")
#   contour_data[[i]]    =  df1
#   power_estimate[[i]]  =  df2
# }
# 
# names(contour_data)  =  names(power_estimate)  =  names(filtered_otu_list)

# 
# 
# #####################Gam Fit Models##############################
# nsamp_vec    =  c(20,seq(50,200, 50))
# gam_model_samp    =  gm  =  list()
# 
# folder_path_mod  =   paste0(path,"gam_model/")
# 
# for(i in 1:7){
#   nam     =   names(filtered_otu_list)[i]    
#   
#   for(j in 1:length(nsamp_vec)){
#     nsam     =   nsamp_vec[j]
#     gm[[j]] =   readRDS(paste0(folder_path_mod,nsam,"_samples/", 
#                                 nam, ".rds"))
# 
#   }
#   
#   names(gm)  = paste0(nsamp_vec, "samples") 
#   gam_model_samp[[i]]    =  gm
# }     
# names(gam_model_samp)  =  names(filtered_otu_list)
# (gam_model_samp[[1]][[1]]$model)
# 
# #####################Gam Fit Models##############################
# notu_vec    = seq(500,1500, 500)
# gam_model_otu    =  gm  =  list()
# 
# folder_path_mod  =   paste0(path,"gam_model/")
# 
# for(i in 1:7){
#   nam     =   names(filtered_otu_list)[i]    
#   
#   for(j in 1:length(notu_vec)){
#     nsam     =   notu_vec[j]
#     gm[[j]] =   readRDS(paste0(folder_path_mod,nsam,"_notu/", 
#                                nam, ".rds"))
#     
#   }
#   
#   names(gm)  = paste0(notu_vec, "otu") 
#   gam_model_otu[[i]]    =  gm
# }     
# names(gam_model_otu)  =  names(filtered_otu_list)
# 

# path_deseq  =   paste0(path,"deseq_results/")
# deseq_sample_size  =  list()
# nsamp_vec  =  seq(10,100,20)
# 
# for(i in 1:7){
#   nam     =   names(filtered_otu_list)[i]
#   dds = foreach(j = 1:length(nsamp_vec),.combine = "rbind") %do% {
#     nsam       =   nsamp_vec[j]
#     dd         =   readRDS(paste0(path_deseq,nsam,"_samples/",nam,".rds"))
#     deseq_list =   do.call(rbind, lapply(dd, function(x) x$deseq_estimate))
#     deseq_list$sample_size = rep(nsam, nrow(deseq_list))
#     deseq_list
#   }
#   deseq_sample_size[[i]]    =  dds
# }
# names(deseq_sample_size)  =  names(filtered_otu_list)

#lapply(countdata_list_obs, dim)

# library(tibble)
# df1 <- enframe(logfoldchange_list) %>% 
#          unnest(value)  %>% 
#         setNames(c("name","logfoldchange")) %>%
#         as.data.frame()
# 
# df2 <- enframe(logmean_list) %>% 
#        unnest(value) %>% 
#        setNames(c("name","logmean")) %>%
#        as.data.frame()
# stopifnot(df1$name == df2$name)
# 
# df  = data.frame(df2,variance=  sqrt(df1$logfoldchange^2))
# 
# oka_col = c(
#   "#E23D28", 
#   "#0000FF",
#   "#E69F00" # Orange
# )
# ggplot(df, aes(logmean, variance)) +
#   geom_point(alpha=0.1)  + 
#   geom_smooth(aes(color="smooth line")) +
#   geom_smooth(method = "lm", formula = y ~poly(x,2), aes(color="quadratic"),lwd=1) +
#   geom_smooth(method = "lm", formula = y ~x, aes(color="linear"),lwd=1) +
#   xlab(TeX("$\\log_2$(mean count)")) + ylab(TeX("|$\\log_2$(fold change)|")) +
#   scale_color_manual(name = "Type",
#                      values = oka_col) +
#   facet_wrap(~name,scale="free")+
#   theme_bw(base_size = 10)+
#   labs(color = " ")#, text = element_text(size = 20))



# scale_location <- function(logmean,logfoldchange){
#   
#   dat <- data.frame(logmean=logmean, 
#                     variance=sqrt(logfoldchange^2))
#   
#   ggplot(dat, aes(logmean, variance)) +
#     geom_point(alpha=0.1)  + 
#     geom_smooth(aes(color="smooth line")) +
#     geom_smooth(method = "lm", formula = y ~poly(x,2), aes(color="quadratic"),lwd=1) +
#     geom_smooth(method = "lm", formula = y ~x, aes(color="linear"),lwd=1) +
#     xlab("logmean") + ylab(expression(sqrt(logfoldchange^2))) +
#     theme_bw(base_size = 10)#, text = element_text(size = 20))
# }




lmean_component_list= readRDS(file = paste0(path,"lmean_component_list.rds"))
lmean_density_plt_list = readRDS(file = paste0(path,"lmean_density_plt_list.rds"))
logmean_param_list=  readRDS(file = paste0(path,"logmean_param_list.rds"))
dispersion_param_list =  readRDS(file = paste0(path,"dispersion_param_list.rds"))
dispersion_confint_list = readRDS(file = paste0(path,"dispersion_confint_list.rds"))
lfc_lm_plot_list  = readRDS(file = paste0(path,"lfc_lm_plot_list.rds"))
optimal_model_list = readRDS(file = paste0(path,"optimal_model_list.rds"))


lfc_plt_list  =  readRDS(file = paste0(path,"lfc_plt_list.rds"))
lfc_aic_list  = readRDS(file = paste0(path,"lfc_aic_list.rds"))
logfoldchange_param_list = readRDS(file = paste0(path,"logfoldchange_param_list.rds"))
countdata_filt_list = readRDS(file = paste0(path,"countdata_filt_list.rds"))
######################################################
source("R_files/Load_Packages.R")
source("R_files/debug_single_example_fun.R")
source("R_files/debug_variance_functions.R")
source("R_files/debug_initial_conditions.R")
######################################################
path     <- "Datasets/"
plt.path <- "Figures/"
######################################################
countdata_list = readRDS(file = "Datasets/data.rds")
metadata_list = readRDS(file = "Datasets/metadata.rds")
logmean_param_list=  readRDS(file = paste0(path,"logmean_param_list.rds"))
dispersion_param_list =  readRDS(file = paste0(path,"dispersion_param_list.rds"))
logfoldchange_param_list = readRDS(file = paste0(path,"logfoldchange_param_list.rds"))
######################################################
n=1; nam    =   names(metadata_list)[n]
countdata =   countdata_list[[n]]
metadata  =   metadata_list[[n]]
logmean_param = logmean_param_list[[n]]
logfoldchange_param = logfoldchange_param_list[[n]]
dispersion_param  =  dispersion_param_list[[n]]
######################################################
v = filter(countdata, metadata,filt_thresh=5)
names(v) 
countdata_obs = v$countdata
#countdat = data.frame(t(v$countdata))
class(countdata_obs)


nsim = 5; nsamp=ncol(countdata_obs)/2; notu=nrow(countdata_obs)
##########################################################
d1 = simulate_fun(logmean_param,logfoldchange_param,
                  dispersion_param,
                  nsamples_per_group=nsamp,
                  notu=notu,nsim=nsim,
                  scale_dispersion=0.3,
                  num_cores = 11,
                  seed = 100)

d1 = d1$countdata_list[1:2]

d2 = simulate_fun(logmean_param,logfoldchange_param,
                  dispersion_param,
                  nsamples_per_group=nsamp,
                  notu=notu,nsim=nsim,
                  scale_dispersion=1,
                  num_cores = 11,
                  seed = 100)

d2 = d2$countdata_list[1:2]
#########################################################



compare_dataset(d1,countdata_obs, data_name= "XYZ")

class(d1)
length(d1)
compare_dataset(countdata_obs,d2)
#########################################################
dlist <- list.append(contdat_list,obs_filt=countdata_filt)

dlist = list.append(contdat_list1,contdat_list,obs_filt=countdata_filt)
dlist <- list.append(contdat_list,obs_filt=countdata_filt)
vars <-(dlist |> purrr::map(~apply(., 1, stats::var)) |>
          purrr::map_dfr(~ tibble(var = .),.id = "type")) 

vars_obs <- vars[vars$type == "obs_filt",]
vars_others <- vars[vars$type != "obs_filt",]

vars_scaled <- vars[vars$type == "scaled",]
vars_unscaled <- vars[vars$type == "unscaled",]
#vars_HMP <- vars[vars$type == "HMP",]

p = ggplot(vars_others, aes(x = log(var), colour = type)) + 
  geom_density(lwd=1.5) +
  geom_density(data=vars_obs, aes(x = log(var)), linetype = "solid",lwd=3) +
  xlab(TeX("$\\log_{10}$(variance of taxa)"))  
p


dlist = list(obs_filt=countdata_obs, scaled=d1,unscaled=d2,HMP=d3)
mean_vec <-(dlist |> purrr::map(~apply(., 1, mean)) |>
              purrr::map_dfr(~ tibble(mean = .),.id = "type")) 

mean_obs <- mean_vec[mean_vec$type == "obs_filt",]
mean_others <- mean_vec[mean_vec$type != "obs_filt",]

p2 = ggplot(mean_others, aes(x = log(mean), colour = type)) + 
  geom_density(lwd=1.5) +
  geom_density(data=mean_obs, aes(x = log(mean)), linetype = "solid",lwd=3) +
  xlab(TeX("$\\log_{10}$(variance of taxa)"))  
p2



# otu_table <- data.frame(t(countdata)) %>%
#   rownames_to_column(var = "Samples")
# 
# combined_otu_table=left_join(metadata,otu_table, by = "Samples")  
# 
# otu_table_control <- combined_otu_table %>%
#   dplyr::filter(Groups == "NT") %>%
#   dplyr::select(-c(Groups)) %>%
#   column_to_rownames("Samples")
# 
# otu_table_treatment <- combined_otu_table %>%
#   dplyr::filter(Groups == "ASD") %>%
#   dplyr::select(-c(Groups)) %>%
#   column_to_rownames("Samples")
# 
# size_cont = rowSums(otu_table_control)
# size_treat = rowSums(otu_table_treatment)


