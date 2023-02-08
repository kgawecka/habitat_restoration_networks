# DISSIMILARITY BETWEEN PATCHES (BETA DIVERSITY)

# this script computes beta diversity between pairs of patches and with the metaweb

library(data.table)
library(dplyr)
library(plyr)
library(betalink)
library(reshape2)
library(parallel)
library(doSNOW)


get_patch_interactions = function(df_out, n_r, n_c, n_patch, M_inc){
  
  # presence data for resources and consumers
  p_res = df_out %>% filter(guild=="resources") %>% select(patch_no, species) %>% 
    dplyr::rename(resource=species) %>%
    mutate(p_res=1)
  p_con = df_out %>% filter(guild=="consumers") %>% select(patch_no, species) %>% 
    dplyr::rename(consumer=species) %>%
    mutate(p_con=1)
  
  # dataframe with all possible interactions
  df = data.frame(expand.grid(resource=1:n_r, consumer=1:n_c, patch_no=1:n_patch),
                  M_inc=rep(as.vector(M_inc), n_patch))
  
  # filter for possible interactions only, then join with presence data
  df_int = filter(df, M_inc==1) %>% 
    left_join(., p_res, by=c("patch_no", "resource")) %>% filter(., p_res==1) %>%
    left_join(., p_con, by=c("patch_no", "consumer")) %>% filter(., p_con==1)
  
  return(df_int)
}

compute_beta_diversity = function(df_int, patches, n_sample, M_inc){
  
  # convert species numbers to factors (must be as factors for betadiversity)
  df_int = df_int %>%
    mutate(resource = as.factor(paste0("r",resource)),
           consumer = as.factor(paste0("c",consumer)))
  
  if(length(patches)>n_sample){
    # sample patches
    p_samp = sample(patches, n_sample)
    df_int_D = df_int %>% filter(patch_no %in% p_samp)
  } else {
    df_int_D = df_int
  }
  
  # create list of local networks
  nets = dlply(df_int_D, .(patch_no), acast, resource~consumer, sum, value.var="M_inc")
  
  # patch 0 = metaweb
  nets$`0` = M_inc
  
  # convert to list of igraph objects
  networks = prepare_networks(nets, directed=TRUE)
  
  # compute beta diversity
  betadiv_current = network_betadiversity(networks)
  
  return(betadiv_current)
}

patch_dissimilarity = function(networkName, networkType, n, dD, RS, n_sample) {
  
  library(data.table)
  library(dplyr)
  library(plyr)
  library(betalink)
  library(reshape2)
  
  # number of patches in grid
  n_patch = n*n
  
  # network incidence matrix
  M_inc = read.table(paste0("../../Data/",networkName,"/M_inc.csv"))
  M_inc = as.matrix(M_inc)
  
  # number of resources and consumers
  n_r = nrow(M_inc)
  n_c = ncol(M_inc)
  rownames(M_inc) = paste0("r",seq(1,n_r))
  colnames(M_inc) = paste0("c",seq(1,n_c))
  
  # DESTRUCTION

  # fractions of habitat loss
  D_all = round(seq(0,1,dD), digits=2)
  
  # initialise dataframe for storing results
  betadivD = data.frame(D=rep(D_all, ((n_sample+1)*n_sample/2)),
                        patch_i=NA,
                        patch_j=NA,
                        beta_S=NA,
                        beta_WN=NA)
  
  for(k in 1:(length(D_all)-1)){
    
    D_current = D_all[k]
    
    # import simulation results
    df_out = fread(paste0("../../Output/",networkName,"/",networkName,"_",networkType,"_D_out_D",floor(D_current*100),".csv"))
    
    if(nrow(df_out)==0 || length(unique(df_out$guild))<2) { break }
    
    # dataframe with patch interactions
    df_int = get_patch_interactions(df_out, n_r, n_c, n_patch, M_inc)
    
    if(nrow(df_int)==0) { next }
    
    patches = unique(df_int$patch_no)
    
    if(length(patches)<=1) { next }
    
    # compute beta diversity
    betadiv_current = compute_beta_diversity(df_int, patches, n_sample, M_inc)
    
    no_results = nrow(betadiv_current)
    
    # store results
    betadivD[betadivD$D==D_current,"patch_i"][1:no_results] = as.integer(levels(betadiv_current$i))[betadiv_current$i]
    betadivD[betadivD$D==D_current,"patch_j"][1:no_results] = as.integer(levels(betadiv_current$j))[betadiv_current$j]
    betadivD[betadivD$D==D_current,"beta_S"][1:no_results] = betadiv_current$S
    betadivD[betadivD$D==D_current,"beta_WN"][1:no_results] = betadiv_current$WN

  } # k loop
  
  # remove empty rows
  betadivD=filter(betadivD, !is.na(patch_i)) %>% mutate(DR=D)

  
  # RESTORATION
  
  # fractions of habitat destruction at which restoration starts
  DR_all = round(seq(dD,(1-dD),dD), digits=2)
  
  # initialise dataframe for storing results
  betadivR = data.frame(expand.grid(DR=DR_all, D=rep(D_all, ((n_sample+1)*n_sample/2))),
                        patch_i=NA,
                        patch_j=NA,
                        beta_S=NA,
                        beta_WN=NA)
  betadivR = filter(betadivR, D<DR)
  
  for(d in 1:length(DR_all)){
    
    DR_current = DR_all[d]
    
    # check if results for DR exist
    if(file.exists(paste0("../../Output/",networkName,"/",networkName,"_",networkType,"_R_",RS,"_DR",floor(DR_current*100),"_out_D0.csv")) == FALSE) { next }
    
    # habitat loss
    D_all = round(seq(DR_current,0,-dD), digits=2)
    
    for(k in 2:length(D_all)){
      
      D_current = D_all[k]
      
      # import simulation results
      df_out = fread(paste0("../../Output/",networkName,"/",networkName,"_",networkType,"_R_",RS,"_DR",floor(DR_current*100),"_out_D",floor(D_current*100),".csv"))
      
      if(nrow(df_out)==0 || length(unique(df_out$guild))<2) { break }
      
      # dataframe with patch interactions
      df_int = get_patch_interactions(df_out, n_r, n_c, n_patch, M_inc)
      
      if(nrow(df_int)==0) { next }
      
      patches = unique(df_int$patch_no)
      
      if(length(patches)<=1) { next }
      
      # compute beta diversity
      betadiv_current = compute_beta_diversity(df_int, patches, n_sample, M_inc)
      
      no_results = nrow(betadiv_current)
      
      # store results
      betadivR[betadivR$DR==DR_current & betadivR$D==D_current,"patch_i"][1:no_results] = as.integer(levels(betadiv_current$i))[betadiv_current$i]
      betadivR[betadivR$DR==DR_current & betadivR$D==D_current,"patch_j"][1:no_results] = as.integer(levels(betadiv_current$j))[betadiv_current$j]
      betadivR[betadivR$DR==DR_current & betadivR$D==D_current,"beta_S"][1:no_results] = betadiv_current$S
      betadivR[betadivR$DR==DR_current & betadivR$D==D_current,"beta_WN"][1:no_results] = betadiv_current$WN
      
    } # k loop
  } # d loop
  
  # remove empty rows
  betadivR=filter(betadivR, !is.na(patch_i))
  
  # combine destruction and restoration dataframes
  betadiv = rbind(betadivR, betadivD)
  
  # write out results
  write.csv(betadiv, paste0("../../Results/",networkName,"/",networkName,"_",networkType,"_R_",RS,"_patch_dissimilarity.csv"), row.names=FALSE)
  
  return(0)
}


# run replicates in parallel
cl = makeCluster(10)  # make a Cluster with 10 cores
#cl <- makeCluster(10, outfile="")  # make a Cluster with the number of cores (Error in unserialize(socklist[[n]]) : error reading from connection)
registerDoSNOW(cl)                    # create Cluster

# Run the simulation for a specified set of values
# This is done in parallel on the created Cluster

#result = foreach(networkName=c("M_SD_025", "M_SD_027", "M_PL_036", "M_SD_008", "M_PL_059", 
#                               "M_SD_014", "M_SD_005", "M_SD_002", "M_PL_033", "M_PL_037")) %dopar% 
#  patch_dissimilarity(networkName, networkType="M", n=100, dD=0.05, RS="diversity", n_sample=100)

#result = foreach(networkName=c("M_PL_025", "M_PL_046", "M_SD_010", "M_SD_012", "M_PL_039", 
#                               "M_PL_006", "M_SD_007", "M_SD_016", "M_PL_051", "M_PL_010")) %dopar% 
#  patch_dissimilarity(networkName, networkType="M", n=100, dD=0.05, RS="diversity", n_sample=100)

#result = foreach(networkName=c("A_HP_015", "A_HP_035", "A_HP_028", "A_HP_005", "A_HP_021", 
#                               "A_HP_007", "A_HP_032", "A_HP_012", "A_HP_008", "A_HP_047")) %dopar% 
#  patch_dissimilarity(networkName, networkType="A", n=100, dD=0.05, RS="nonrandom", n_sample=100)

#result = foreach(networkName=c("A_HP_002", "A_HP_033", "A_HP_042", "A_HP_046", "A_HP_025", 
#                               "A_HP_050", "A_PH_007", "A_PH_004", "A_PH_005", "A_PH_006")) %dopar% 
#  patch_dissimilarity(networkName, networkType="A", n=100, dD=0.05, RS="nonrandom", n_sample=100)

stopCluster(cl)               # when finished stop the Cluster
