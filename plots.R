source("Load_Packages.R")
source("simulate_fun.R")
source("utils.R")
source("fitting_fun.R")
source("read_estimate_results.R")
########################Test####################################################
path  =  "data/"
countdata_list_obs  =  readRDS(file = paste0(path,"autism_data.rds"))
metadata_list_obs   =  readRDS(file = paste0(path,"autism_metadata.rds"))
################################################################################
#####Show examples of log(fold change) vrs log(mean abundance) plot
c = c(2,5,3)
dd = foreach(i = c,.combine = "rbind") %do% {
  nam             =   names(logmean_list)[i]
  logmean         =   logmean_list[[i]]
  logfoldchange   =   logfoldchange_list[[i]]
         dd       =   data.frame(logmean, logfoldchange)
  dd$datatype     =   rep(nam, nrow(dd))
  dd
}

ggplot(dd, aes(logmean, logfoldchange)) +
  geom_point(alpha=0.1) +
  geom_smooth() +
  xlab(TeX("$\\log_2$(mean abundance)")) +
  ylab(TeX("$\\log_2$(fold change)")) +
  facet_wrap(~datatype, scale="free") 
#########################################################################
##Compare simulations from HMP, metaSPARSim and our method
res  = readRDS(paste0(path,"HMP_metaSPARSim_sim_compare.rds")) 

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
  countdata_sim     =    list(HMP = HMP, metaSPARSim = metaSPARSim,
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
###############################################################################
nsamp      =   100
path1      =   paste0(path,"contour_data/",nsamp,"_samples/")
pltt   =   list()
for(i in 1:7){
  nam            =   names(dispersion_param_list)[i]
  combined_data  =   readRDS(paste0(path1,"combined_data_",nam,"_.rds"))
  power_estimate =   readRDS(paste0(path1,"power_estimate_",nam,"_.rds"))
  
  combined_data$pvalue_reject <- factor(combined_data$pval_reject)
  combined_data$pvalue_reject <- factor(ifelse(combined_data$pvalue_reject == 0,
                                               "No", "Yes"))
  
  cont_breaks =  seq(0,1,0.1)
  pltt[[i]] <- (ggplot(combined_data)
                    + aes(lmean_abund, abs_lfc)
                    + rasterise(geom_point(aes(color = pvalue_reject), alpha = 0.5))
                    + xlab(TeX("$\\log_2$(mean counts)"))
                    + ylab(TeX("|$\\log_2$(fold change)|"))
                    + scale_colour_manual(values = c("black", "red"))
                    + labs(color="reject null hypothesis")
                    + geom_contour(data = power_estimate,
                                   aes(z=power),lwd=1,
                                   breaks = cont_breaks)
                    + geom_label_contour(data = power_estimate,
                                         aes(z= power,label = sprintf("%.3f", after_stat(level))),
                                         breaks = cont_breaks
                    )
  )
  
}

(pltt[[1]] | pltt[[2]] | pltt[[3]]) /
(pltt[[4]] | pltt[[5]] | pltt[[6]]) /
(pltt[[7]] | plot_spacer() | plot_spacer()) +  plot_layout(guides = "collect")
###############################################################################
## Relationship between statistical power, sample size and log fold change 
## (log mean abundance = 5)

nsamples =  seq(30,190,20)
folder_path  =   paste0(path,"gam_model/")

##plot quantiles and average power
#' simulate data and predict 
nsamples =  seq(30,190,20)
notu = 1000; disp_scale = 0.3; nsim= 50; nsamp = 50
folder_path  =   paste0(path,"gam_model/")

otu_data_list  = list()
for(i in 1:7){
  dispersion_param      =   dispersion_param_list[[i]]
  logmean_param         =   logmean_param_list[[i]]
  logfoldchange_param   =   logfoldchange_param_list[[i]]
  otu_data_list[[i]]   =  countdata_sim_fun(logmean_param, logfoldchange_param, 
                                            dispersion_param,
                                            nsamp_per_group=nsamp,
                                            disp_scale = disp_scale,
                                            ncont = NULL,ntreat = NULL,
                                            notu, nsim, 
                                            maxlfc_iter = 10000,
                                            seed = 100)
}

##########################################################################
pp = foreach(i =  1:7, .combine = "rbind") %:%
  foreach(j =  1:length(nsamples),.combine = "rbind") %do%{
    
    nam        =   names(metadata_list_obs)[i]
    otu_data   =   otu_data_list[[i]]
    
    nsam       =   nsamples[j]
      mod      =   readRDS(paste0(folder_path,nsam,"_samples/",
                              nam, ".rds"))
  
     newdata   =   data.frame(lmean_abund  =  unlist(otu_data$logmean_list),
                              abs_lfc      =  unlist(otu_data$logfoldchange_list)) 

    power      =   predict(mod,newdata,type="response")
  quantiles    =   quantile(power, probs = c(0.5, 0.6, 0.7, 0.8))
  
    ddf        =   data.frame(sample_size = nsam, average_power = mean(power))
    ddf        =   cbind(ddf, as.data.frame(t(quantiles))) 
    
  # Rename the quantile columns
  colnames(ddf)[3:6] <- c("quantile_50%", "quantile_60%", "quantile_70%", "quantile_80%")
  ddf$name =   nam
  ddf
}

#oka_col = c("#556B2F", "#E23D28", "#0000FF","#000000")

oka_col = c(
  "#556B2F", 
  "#E23D28", 
  "#0000FF",
  "#E69F00", 
  "#000000"
)


sub = c("PRJNA168470","PRJNA355023","PRJNA644763","PRJNA589343")
pp_sub = pp[pp$name %in% sub,]

ggplot(pp_sub, aes(sample_size, average_power)) +
  geom_point(aes(shape = "average power across all taxa"), size = 3, color = "black") +  # Use shape for points
  geom_line(aes(y = `quantile_50%`, color = "50th Quantile", group = 1), lwd = 1.5) +
  geom_line(aes(y = `quantile_60%`, color = "60th Quantile", group = 1), lwd = 1.5) +
  geom_line(aes(y = `quantile_70%`, color = "70th Quantile", group = 1), lwd = 1.5) +
  geom_line(aes(y = `quantile_80%`, color = "80th Quantile", group = 1), lwd = 1.5) +
  scale_color_manual(name = "Quantile",  # Color legend for lines
                     values = c("50th Quantile" = oka_col[1], 
                                "60th Quantile" = oka_col[2], 
                                "70th Quantile" = oka_col[3], 
                                "80th Quantile" = oka_col[4])) +
  scale_shape_manual(name = "Average Power",  # Separate shape legend for points
                     values = c("average power across all taxa" = 16)) +  # Choose a shape for points
  custom_theme(16) + 
  labs(x = "Sample Size", y = "Average Power") +
  facet_wrap(~name, scales = "free") +
  guides(color = guide_legend(order = 1, reverse = TRUE),  # Reverse color legend
         shape = guide_legend(order = 2))  

###############################################################################
##' Expected number of significant taxa

oka_col_cvd = c(
  "#000000", # Black
  "#E69F00", # Orange
  "#56B4E9", # Sky Blue
  "#009E73", # Bluish Green
  "#F0E442", # Yellow
  "#0072B2", # Blue
  "#D55E00"  # Vermilion
)

pp$expected  = notu * pp$average_power

ggplot(pp, aes(sample_size,expected, group=name, color=name)) +
  geom_point(size=3) +
  geom_line(lwd=1.2)  +
  scale_color_manual(name = "data accession no.", 
                     values= oka_col_cvd)+
  custom_theme(15) + labs(x="sample size", 
                          y="expected number of significant taxa") 

######################################################################
## Relationship between sample size, power, effect size and log fold change
newdata    =   data.frame(lmean_abund  = c(5,5,5),
                          abs_lfc      =  c(2,3,4)) 

pp1 = foreach(i =  1:7, .combine = "rbind") %:%
  foreach(j =  1:length(nsamples),.combine = "rbind") %do%{
    
    nam        =   names(metadata_list_obs)[i]
    nsam       =   nsamples[j]
    mod        =   readRDS(paste0(folder_path,nsam,"_samples/",
                                nam, ".rds"))
    
    power      =   predict(mod,newdata,type="response")

    ddf        =   data.frame(sample_size = rep(nsam,nrow(newdata)), 
                              power = power)
    ddf        =   cbind(newdata,ddf)
  ddf$abs_lfc  =   factor(ddf$abs_lfc, levels = unique(ddf$abs_lfc))
    ddf$name   =   rep(nam,nrow(newdata))
    ddf
  }


sub = c("PRJNA168470","PRJNA687773","PRJNA644763","PRJNA453621")
pp1_sub = pp1[pp1$name %in% sub,]

oka_col_cvd = c(
  "#E69F00", # Orange
  "#56B4E9", # Sky Blue
  "#009E73" # Bluish Green
)

ggplot(pp1_sub, aes(sample_size,power, group=abs_lfc, color=abs_lfc)) +
  geom_point(size=3) +
  geom_line(lwd=1.2)  +
 scale_color_manual(name = TeX("$\\log_{2}$(fold change)"), 
                     values= oka_col_cvd)+
  custom_theme(18) + labs(x="sample size", 
                          y="power") +
  facet_wrap(~name) +
  guides(color = guide_legend(reverse = TRUE), 
         linetype = guide_legend(reverse = TRUE))
################################################################################
# xlab(TeX("$\\log_2$(mean abundance)"))
