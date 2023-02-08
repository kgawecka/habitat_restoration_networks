# COMBINE ALL NETWORK RESULTS

# load packages
library(dplyr)
library(tidyr)
library(data.table)
library(igraph)

# set working directory to "habitat_restoration_networks"
setwd("~/habitat_restoration_networks")


networks = c("M_SD_025", "M_SD_027", "M_PL_036", "M_SD_008", "M_PL_059",
             "M_SD_014", "M_SD_005", "M_SD_002", "M_PL_033", "M_PL_037",
             "M_PL_025", "M_PL_046", "M_SD_010", "M_SD_012", "M_PL_039",
             "M_PL_006", "M_SD_007", "M_SD_016", "M_PL_051", "M_PL_010")


# SPECIES DEGREE ----

sp_deg = data.frame(network=factor(),
                    guild=factor(),
                    species=integer(),
                    degree=integer(),
                    degree_fract=double())

for(network in networks){
  
  M_inc = read.table(paste0("Data/",network,"/M_inc.csv"), quote="\"", comment.char="")
  n_r = nrow(M_inc)
  n_c = ncol(M_inc)
  
  deg = rbind(data.frame(network=as.factor(network),
                         guild="resources",
                         species=1:n_r,
                         degree=rowSums(M_inc),
                         degree_fract=rowSums(M_inc)/n_c),
              data.frame(network=as.factor(network),
                         guild="consumers",
                         species=1:n_c,
                         degree=colSums(M_inc),
                         degree_fract=colSums(M_inc)/n_r))
  
  sp_deg = rbind(sp_deg, deg)
  
}

# write out results
write.csv(sp_deg, paste0("Results/networks_combined_M/species_degree.csv", sep=""), row.names=FALSE)


# INTERACTION DEGREE ----

species_degree = fread(paste0("Results/networks_combined_M/species_degree.csv"))

int_deg = data.frame(network=factor(),
                     interaction=integer(),
                     resource=integer(),
                     consumer=integer(),
                     degree=integer(),
                     degree_fract=double())

for(net in networks){
  
  intID = fread(paste0("Data/",net,"/interactionID.csv"))
  
  deg = intID %>%
    left_join(., filter(species_degree, network==net, guild=="resources"), by=c("resource"="species")) %>%
    left_join(., filter(species_degree, network==net, guild=="consumers"), by=c("consumer"="species")) %>%
    mutate(degree=degree.x+degree.y, degree_fract=(degree_fract.x+degree_fract.y)/2) %>%
    select(network.x, interaction, resource, consumer, degree, degree_fract) %>%
    dplyr::rename(network=network.x)
  
  int_deg = rbind(int_deg, deg)
  
}

# write out results
write.csv(int_deg, paste0("Results/networks_combined_M/interaction_degree.csv", sep=""), row.names=FALSE)


# ABUNDANCE INTERACTIONS ----

abundanceD = data.frame(D=double(),
                        interaction=factor(),
                        abundance=double(),
                        network=factor())

for(network in networks){
  
  ab = fread(paste0("Results/",network,"/",network,"_M_D_abundance_interactions.csv")) %>%
    mutate(network=as.factor(network))
  
  abundanceD = rbind(abundanceD, ab)
  
}

# write out results
write.csv(abundanceD, paste0("Results/networks_combined_M/abundance_interactionsD.csv", sep=""), row.names=FALSE)


abundanceR = data.frame(DR=double(),
                        D=double(),
                        interaction=factor(),
                        abundance=double(),
                        network=factor(),
                        strategy=factor())

for(network in networks){
  
  ab = fread(paste0("Results/",network,"/",network,"_M_R_random_abundance_interactions.csv")) %>%
    mutate(network=as.factor(network), strategy=as.factor("random"))
  
  abundanceR = rbind(abundanceR, ab)
  
  ab = fread(paste0("Results/",network,"/",network,"_M_R_nonrandom_abundance_interactions.csv")) %>%
    mutate(network=as.factor(network), strategy=as.factor("nonrandom"))
  
  abundanceR = rbind(abundanceR, ab)
  
  ab = fread(paste0("Results/",network,"/",network,"_M_R_diversity_abundance_interactions.csv")) %>%
    mutate(network=as.factor(network), strategy=as.factor("diversity"))
  
  abundanceR = rbind(abundanceR, ab)
  
  ab = fread(paste0("Results/",network,"/",network,"_M_R_generalist_abundance_interactions.csv")) %>%
    mutate(network=as.factor(network), strategy=as.factor("generalist"))
  
  abundanceR = rbind(abundanceR, ab)
  
  ab = fread(paste0("Results/",network,"/",network,"_M_R_specialist_abundance_interactions.csv")) %>%
    mutate(network=as.factor(network), strategy=as.factor("specialist"))
  
  abundanceR = rbind(abundanceR, ab)
  
}

# write out results
write.csv(abundanceR, paste0("Results/networks_combined_M/abundance_interactionsR.csv", sep=""), row.names=FALSE)


# EXTINCTION THRESHOLDS ----

abundanceD = fread(paste0("Results/networks_combined_M/abundance_interactionsD.csv"))

D_ext = abundanceD %>%
  filter(abundance>0) %>%
  dplyr::group_by(network, interaction) %>%
  dplyr::summarise(De=max(D)+0.05) %>%
  ungroup()

# write out results
write.csv(D_ext, paste0("Results/networks_combined_M/extinction_threshold_interactions.csv", sep=""), row.names=FALSE)


# RESTORATION EFFICIENCY INTERACTIONS ----

RE = data.frame(DR=double(),
                interaction=factor(),
                RE=double(),
                network=factor(),
                strategy=factor())

for(network in networks){
  
  RE_out = fread(paste0("Results/",network,"/",network,"_M_random_restoration_efficiency_interactions.csv")) %>%
    mutate(network=as.factor(network), strategy=as.factor("random"))
  
  RE = rbind(RE, RE_out)
  
  RE_out = fread(paste0("Results/",network,"/",network,"_M_nonrandom_restoration_efficiency_interactions.csv")) %>%
    mutate(network=as.factor(network), strategy=as.factor("nonrandom"))
  
  RE = rbind(RE, RE_out)
  
  RE_out = fread(paste0("Results/",network,"/",network,"_M_diversity_restoration_efficiency_interactions.csv")) %>%
    mutate(network=as.factor(network), strategy=as.factor("diversity"))
  
  RE = rbind(RE, RE_out)
  
  RE_out = fread(paste0("Results/",network,"/",network,"_M_generalist_restoration_efficiency_interactions.csv")) %>%
    mutate(network=as.factor(network), strategy=as.factor("generalist"))
  
  RE = rbind(RE, RE_out)
  
  RE_out = fread(paste0("Results/",network,"/",network,"_M_specialist_restoration_efficiency_interactions.csv")) %>%
    mutate(network=as.factor(network), strategy=as.factor("specialist"))
  
  RE = rbind(RE, RE_out)
  
}

# write out results
write.csv(RE, paste0("Results/networks_combined_M/restoration_efficiency_interactions.csv", sep=""), row.names=FALSE)


# PATCH NETWORK SIZE ----

patch_size_summary = data.frame(DR=double(),
                                D=double(),
                                n_sp_mean=double(),
                                n_sp_median=double(),
                                n_sp_sd=double(),
                                n_sp_IQR=double(),
                                fract_sp_mean=double(),
                                fract_sp_median=double(),
                                fract_sp_sd=double(),
                                fract_sp_IQR=double(),
                                n_int_mean=double(),
                                n_int_median=double(),
                                n_int_sd=double(),
                                n_int_IQR=double(),
                                fract_int_mean=double(),
                                fract_int_median=double(),
                                fract_int_sd=double(),
                                fract_int_IQR=double(),
                                network=factor(),
                                strategy=factor())

for(network in networks){
  
  patch_size = fread(paste0("Results/",network,"/",network,"_M_R_random_patch_size.csv"))

  df_summary = patch_size %>% 
    dplyr::group_by(DR, D) %>%
    dplyr::summarise(n_sp_mean=mean(n_sp, na.rm=TRUE), n_sp_median=median(n_sp, na.rm=TRUE),
                     n_sp_sd=sd(n_sp, na.rm=TRUE), n_sp_IQR=IQR(n_sp, na.rm=TRUE),
                     fract_sp_mean=mean(fract_sp, na.rm=TRUE), fract_sp_median=median(fract_sp, na.rm=TRUE),
                     fract_sp_sd=sd(fract_sp, na.rm=TRUE), fract_sp_IQR=IQR(fract_sp, na.rm=TRUE),
                     n_int_mean=mean(n_int, na.rm=TRUE), n_int_median=median(n_int, na.rm=TRUE),
                     n_int_sd=sd(n_int, na.rm=TRUE), n_int_IQR=IQR(n_int, na.rm=TRUE),
                     fract_int_mean=mean(fract_int, na.rm=TRUE), fract_int_median=median(fract_int, na.rm=TRUE),
                     fract_int_sd=sd(fract_int, na.rm=TRUE), fract_int_IQR=IQR(fract_int, na.rm=TRUE)) %>%
    ungroup() %>%
    mutate(network=as.factor(network), strategy=as.factor("random"))
  
  patch_size_summary = rbind(patch_size_summary, df_summary)
  
  
  patch_size = fread(paste0("Results/",network,"/",network,"_M_R_nonrandom_patch_size.csv"))
  
  df_summary = patch_size %>% 
    dplyr::group_by(DR, D) %>%
    dplyr::summarise(n_sp_mean=mean(n_sp, na.rm=TRUE), n_sp_median=median(n_sp, na.rm=TRUE),
                     n_sp_sd=sd(n_sp, na.rm=TRUE), n_sp_IQR=IQR(n_sp, na.rm=TRUE),
                     fract_sp_mean=mean(fract_sp, na.rm=TRUE), fract_sp_median=median(fract_sp, na.rm=TRUE),
                     fract_sp_sd=sd(fract_sp, na.rm=TRUE), fract_sp_IQR=IQR(fract_sp, na.rm=TRUE),
                     n_int_mean=mean(n_int, na.rm=TRUE), n_int_median=median(n_int, na.rm=TRUE),
                     n_int_sd=sd(n_int, na.rm=TRUE), n_int_IQR=IQR(n_int, na.rm=TRUE),
                     fract_int_mean=mean(fract_int, na.rm=TRUE), fract_int_median=median(fract_int, na.rm=TRUE),
                     fract_int_sd=sd(fract_int, na.rm=TRUE), fract_int_IQR=IQR(fract_int, na.rm=TRUE)) %>%
    ungroup() %>%
    mutate(network=as.factor(network), strategy=as.factor("nonrandom"))
  
  patch_size_summary = rbind(patch_size_summary, df_summary)
  
  
  patch_size = fread(paste0("Results/",network,"/",network,"_M_R_diversity_patch_size.csv"))
  
  df_summary = patch_size %>% 
    dplyr::group_by(DR, D) %>%
    dplyr::summarise(n_sp_mean=mean(n_sp, na.rm=TRUE), n_sp_median=median(n_sp, na.rm=TRUE),
                     n_sp_sd=sd(n_sp, na.rm=TRUE), n_sp_IQR=IQR(n_sp, na.rm=TRUE),
                     fract_sp_mean=mean(fract_sp, na.rm=TRUE), fract_sp_median=median(fract_sp, na.rm=TRUE),
                     fract_sp_sd=sd(fract_sp, na.rm=TRUE), fract_sp_IQR=IQR(fract_sp, na.rm=TRUE),
                     n_int_mean=mean(n_int, na.rm=TRUE), n_int_median=median(n_int, na.rm=TRUE),
                     n_int_sd=sd(n_int, na.rm=TRUE), n_int_IQR=IQR(n_int, na.rm=TRUE),
                     fract_int_mean=mean(fract_int, na.rm=TRUE), fract_int_median=median(fract_int, na.rm=TRUE),
                     fract_int_sd=sd(fract_int, na.rm=TRUE), fract_int_IQR=IQR(fract_int, na.rm=TRUE)) %>%
    ungroup() %>%
    mutate(network=as.factor(network), strategy=as.factor("diversity"))
  
  patch_size_summary = rbind(patch_size_summary, df_summary)
  
  
  patch_size = fread(paste0("Results/",network,"/",network,"_M_R_specialist_patch_size.csv"))

  df_summary = patch_size %>% 
    dplyr::group_by(DR, D) %>%
    dplyr::summarise(n_sp_mean=mean(n_sp, na.rm=TRUE), n_sp_median=median(n_sp, na.rm=TRUE),
                     n_sp_sd=sd(n_sp, na.rm=TRUE), n_sp_IQR=IQR(n_sp, na.rm=TRUE),
                     fract_sp_mean=mean(fract_sp, na.rm=TRUE), fract_sp_median=median(fract_sp, na.rm=TRUE),
                     fract_sp_sd=sd(fract_sp, na.rm=TRUE), fract_sp_IQR=IQR(fract_sp, na.rm=TRUE),
                     n_int_mean=mean(n_int, na.rm=TRUE), n_int_median=median(n_int, na.rm=TRUE),
                     n_int_sd=sd(n_int, na.rm=TRUE), n_int_IQR=IQR(n_int, na.rm=TRUE),
                     fract_int_mean=mean(fract_int, na.rm=TRUE), fract_int_median=median(fract_int, na.rm=TRUE),
                     fract_int_sd=sd(fract_int, na.rm=TRUE), fract_int_IQR=IQR(fract_int, na.rm=TRUE)) %>%
    ungroup() %>%
    mutate(network=as.factor(network), strategy=as.factor("specialist"))
  
  patch_size_summary = rbind(patch_size_summary, df_summary)
  
  
  patch_size = fread(paste0("Results/",network,"/",network,"_M_R_generalist_patch_size.csv"))

  df_summary = patch_size %>% 
    dplyr::group_by(DR, D) %>%
    dplyr::summarise(n_sp_mean=mean(n_sp, na.rm=TRUE), n_sp_median=median(n_sp, na.rm=TRUE),
                     n_sp_sd=sd(n_sp, na.rm=TRUE), n_sp_IQR=IQR(n_sp, na.rm=TRUE),
                     fract_sp_mean=mean(fract_sp, na.rm=TRUE), fract_sp_median=median(fract_sp, na.rm=TRUE),
                     fract_sp_sd=sd(fract_sp, na.rm=TRUE), fract_sp_IQR=IQR(fract_sp, na.rm=TRUE),
                     n_int_mean=mean(n_int, na.rm=TRUE), n_int_median=median(n_int, na.rm=TRUE),
                     n_int_sd=sd(n_int, na.rm=TRUE), n_int_IQR=IQR(n_int, na.rm=TRUE),
                     fract_int_mean=mean(fract_int, na.rm=TRUE), fract_int_median=median(fract_int, na.rm=TRUE),
                     fract_int_sd=sd(fract_int, na.rm=TRUE), fract_int_IQR=IQR(fract_int, na.rm=TRUE)) %>%
    ungroup() %>%
    mutate(network=as.factor(network), strategy=as.factor("generalist"))
  
  patch_size_summary = rbind(patch_size_summary, df_summary)
  
}

# write out results
write.csv(patch_size_summary, paste0("Results/networks_combined_M/patch_size_summary.csv", sep=""), row.names=FALSE)


# PATCH NETWORK DISSIMILARITY ----

beta_summary_patch = data.frame(DR=double(),
                                D=double(),
                                beta_S_mean=double(),
                                beta_S_median=double(),
                                beta_S_sd=double(),
                                beta_S_IQR=double(),
                                beta_WN_mean=double(),
                                beta_WN_median=double(),
                                beta_WN_sd=double(),
                                beta_WN_IQR=double(),
                                network=factor(),
                                strategy=factor())

beta_summary_meta = data.frame(DR=double(),
                               D=double(),
                               beta_S_mean=double(),
                               beta_S_median=double(),
                               beta_S_sd=double(),
                               beta_S_IQR=double(),
                               beta_WN_mean=double(),
                               beta_WN_median=double(),
                               beta_WN_sd=double(),
                               beta_WN_IQR=double(),
                               network=factor(),
                               strategy=factor())

for(network in networks){
  
  betadiv = fread(paste0("Results/",network,"/",network,"_M_R_random_patch_dissimilarity.csv"))
  
  df_summary = betadiv %>% 
    filter(patch_j>0) %>%
    dplyr::group_by(DR, D) %>%
    dplyr::summarise(beta_S_mean=mean(beta_S, na.rm=TRUE), beta_S_median=median(beta_S, na.rm=TRUE),
                     beta_S_sd=sd(beta_S, na.rm=TRUE), beta_S_IQR=IQR(beta_S, na.rm=TRUE),
                     beta_WN_mean=mean(beta_WN, na.rm=TRUE), beta_WN_median=median(beta_WN, na.rm=TRUE),
                     beta_WN_sd=sd(beta_WN, na.rm=TRUE), beta_WN_IQR=IQR(beta_WN, na.rm=TRUE)) %>%
    ungroup() %>%
    mutate(network=as.factor(network), strategy=as.factor("random"))
  
  beta_summary_patch = rbind(beta_summary_patch, df_summary)
  
  df_summary = betadiv %>% 
    filter(patch_j==0) %>%
    dplyr::group_by(DR, D) %>%
    dplyr::summarise(beta_S_mean=mean(beta_S, na.rm=TRUE), beta_S_median=median(beta_S, na.rm=TRUE),
                     beta_S_sd=sd(beta_S, na.rm=TRUE), beta_S_IQR=IQR(beta_S, na.rm=TRUE),
                     beta_WN_mean=mean(beta_WN, na.rm=TRUE), beta_WN_median=median(beta_WN, na.rm=TRUE),
                     beta_WN_sd=sd(beta_WN, na.rm=TRUE), beta_WN_IQR=IQR(beta_WN, na.rm=TRUE)) %>%
    ungroup() %>%
    mutate(network=as.factor(network), strategy=as.factor("random"))
  
  beta_summary_meta = rbind(beta_summary_meta, df_summary)
  
  
  betadiv = fread(paste0("Results/",network,"/",network,"_M_R_nonrandom_patch_dissimilarity.csv"))
  
  df_summary = betadiv %>% 
    filter(patch_j>0) %>%
    dplyr::group_by(DR, D) %>%
    dplyr::summarise(beta_S_mean=mean(beta_S, na.rm=TRUE), beta_S_median=median(beta_S, na.rm=TRUE),
                     beta_S_sd=sd(beta_S, na.rm=TRUE), beta_S_IQR=IQR(beta_S, na.rm=TRUE),
                     beta_WN_mean=mean(beta_WN, na.rm=TRUE), beta_WN_median=median(beta_WN, na.rm=TRUE),
                     beta_WN_sd=sd(beta_WN, na.rm=TRUE), beta_WN_IQR=IQR(beta_WN, na.rm=TRUE)) %>%
    ungroup() %>%
    mutate(network=as.factor(network), strategy=as.factor("nonrandom"))
  
  beta_summary_patch = rbind(beta_summary_patch, df_summary)
  
  df_summary = betadiv %>% 
    filter(patch_j==0) %>%
    dplyr::group_by(DR, D) %>%
    dplyr::summarise(beta_S_mean=mean(beta_S, na.rm=TRUE), beta_S_median=median(beta_S, na.rm=TRUE),
                     beta_S_sd=sd(beta_S, na.rm=TRUE), beta_S_IQR=IQR(beta_S, na.rm=TRUE),
                     beta_WN_mean=mean(beta_WN, na.rm=TRUE), beta_WN_median=median(beta_WN, na.rm=TRUE),
                     beta_WN_sd=sd(beta_WN, na.rm=TRUE), beta_WN_IQR=IQR(beta_WN, na.rm=TRUE)) %>%
    ungroup() %>%
    mutate(network=as.factor(network), strategy=as.factor("nonrandom"))
  
  beta_summary_meta = rbind(beta_summary_meta, df_summary)
  
  
  betadiv = fread(paste0("Results/",network,"/",network,"_M_R_diversity_patch_dissimilarity.csv"))
  
  df_summary = betadiv %>% 
    filter(patch_j>0) %>%
    dplyr::group_by(DR, D) %>%
    dplyr::summarise(beta_S_mean=mean(beta_S, na.rm=TRUE), beta_S_median=median(beta_S, na.rm=TRUE),
                     beta_S_sd=sd(beta_S, na.rm=TRUE), beta_S_IQR=IQR(beta_S, na.rm=TRUE),
                     beta_WN_mean=mean(beta_WN, na.rm=TRUE), beta_WN_median=median(beta_WN, na.rm=TRUE),
                     beta_WN_sd=sd(beta_WN, na.rm=TRUE), beta_WN_IQR=IQR(beta_WN, na.rm=TRUE)) %>%
    ungroup() %>%
    mutate(network=as.factor(network), strategy=as.factor("diversity"))
  
  beta_summary_patch = rbind(beta_summary_patch, df_summary)
  
  df_summary = betadiv %>% 
    filter(patch_j==0) %>%
    dplyr::group_by(DR, D) %>%
    dplyr::summarise(beta_S_mean=mean(beta_S, na.rm=TRUE), beta_S_median=median(beta_S, na.rm=TRUE),
                     beta_S_sd=sd(beta_S, na.rm=TRUE), beta_S_IQR=IQR(beta_S, na.rm=TRUE),
                     beta_WN_mean=mean(beta_WN, na.rm=TRUE), beta_WN_median=median(beta_WN, na.rm=TRUE),
                     beta_WN_sd=sd(beta_WN, na.rm=TRUE), beta_WN_IQR=IQR(beta_WN, na.rm=TRUE)) %>%
    ungroup() %>%
    mutate(network=as.factor(network), strategy=as.factor("diversity"))
  
  beta_summary_meta = rbind(beta_summary_meta, df_summary)
  
}

# write out results
write.csv(beta_summary_patch, paste0("Results/networks_combined_M/dissimilarity_patch_summary.csv", sep=""), row.names=FALSE)
write.csv(beta_summary_meta, paste0("Results/networks_combined_M/dissimilarity_meta_summary.csv", sep=""), row.names=FALSE)


# RESTORATION EFFICIENCY NETWORKS ----

RE = data.frame(DR=double(),
                RE_size_sp=double(),
                RE_size_int=double(),
                RE_beta_patch_sp=double(),
                RE_beta_patch_int=double(),
                RE_beta_meta_sp=double(),
                RE_beta_meta_int=double(),
                network=factor(),
                strategy=factor())

for(network in networks){
  
  RE_out = fread(paste0("Results/",network,"/",network,"_M_random_restoration_efficiency_networks.csv")) %>%
    mutate(network=as.factor(network), strategy=as.factor("random"))
  
  RE = rbind(RE, RE_out)
  
  RE_out = fread(paste0("Results/",network,"/",network,"_M_nonrandom_restoration_efficiency_networks.csv")) %>%
    mutate(network=as.factor(network), strategy=as.factor("nonrandom"))
  
  RE = rbind(RE, RE_out)
  
  RE_out = fread(paste0("Results/",network,"/",network,"_M_diversity_restoration_efficiency_networks.csv")) %>%
    mutate(network=as.factor(network), strategy=as.factor("diversity"))
  
  RE = rbind(RE, RE_out)
  
  RE_out = fread(paste0("Results/",network,"/",network,"_M_specialist_restoration_efficiency_networks.csv")) %>%
    mutate(network=as.factor(network), strategy=as.factor("specialist"))
  
  RE = rbind(RE, RE_out)
  
  RE_out = fread(paste0("Results/",network,"/",network,"_M_generalist_restoration_efficiency_networks.csv")) %>%
    mutate(network=as.factor(network), strategy=as.factor("generalist"))
  
  RE = rbind(RE, RE_out)

}

# write out results
write.csv(RE, paste0("Results/networks_combined_M/restoration_efficiency_networks.csv", sep=""), row.names=FALSE)