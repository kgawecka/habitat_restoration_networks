# METANETWORK STRUCTURE

library(data.table)
library(dplyr)
library(igraph)
library(FactoMineR)

setwd("~/habitat_restoration_networks")


networks = c("M_SD_005", "M_SD_008", "M_SD_010", "M_SD_012", "M_SD_025",
             "M_SD_002", "M_SD_007", "M_SD_014", "M_SD_016", "M_SD_027",
             "M_PL_006", "M_PL_010", "M_PL_036", "M_PL_037", "M_PL_059",
             "M_PL_025", "M_PL_033", "M_PL_039", "M_PL_046", "M_PL_051",
             "A_HP_015", "A_HP_035", "A_HP_028", "A_HP_005", "A_HP_021", 
             "A_HP_007", "A_HP_032", "A_HP_012", "A_HP_008", "A_HP_047", 
             "A_HP_002", "A_HP_033", "A_HP_042", "A_HP_046", "A_HP_025",
             "A_HP_050", "A_PH_007", "A_PH_004", "A_PH_005", "A_PH_006")

compute_nestedness = function(B){
  
  # Get number of rows and columns
  nrows <- nrow(B)
  ncols <- ncol(B)
  
  # Compute nestedness of rows
  nestedness_rows <- 0
  for(i in 1:(nrows-1)){
    for(j in (i+1): nrows){
      
      c_ij <- sum(B[i,] * B[j,])      # Number of interactions shared by i and j
      k_i <- sum(B[i,])               # Degree of node i
      k_j <- sum(B[j,])               # Degree of node j
      
      if (k_i == 0 || k_j==0) {next}  # Handle case if a node is disconnected
      
      o_ij <- c_ij / min(k_i, k_j)    # Overlap between i and j
      
      nestedness_rows <- nestedness_rows + o_ij
    }
  }
  
  # Compute nestedness of columns
  nestedness_cols <- 0
  for(i in 1: (ncols-1)){
    for(j in (i+1): ncols){
      
      c_ij <- sum(B[,i] * B[,j])      # Number of interactions shared by i and j
      k_i <- sum(B[,i])               # Degree of node i
      k_j <- sum(B[,j])               # Degree of node j
      if (k_i == 0 || k_j==0) {next}  # Handle case if a node is disconnected.
      
      o_ij <- c_ij / min(k_i, k_j)    # Overlap between i and j
      
      nestedness_cols <- nestedness_cols + o_ij         
    }
  }
  
  # Compute nestedness of the network
  nestedness <- (nestedness_rows + nestedness_cols) / ((nrows * (nrows - 1) / 2) + (ncols * (ncols - 1) / 2))
  
  return(nestedness)
}

cell_model = function(d,t_max){
  
  rows = nrow(d)
  columns = ncol(d)
  
  t=1
  
  null_nest = rep(NA,t_max)
  null_mod = rep(NA,t_max)
  
  while (t <= t_max){
    
    PR = matrix(0, rows, 1)
    PC = matrix(0, columns, 1)
    B = matrix(0, rows, columns)
    
    for (i in 1:rows){
      number_ones=0
      for (j in 1:columns){
        if(d[i,j] == 1){
          number_ones=number_ones+1
        }
      }
      PR[i] = number_ones/columns
    }
    
    for (j in 1:columns){
      number_ones=0
      for (i in 1:rows){
        if(d[i,j] == 1){
          number_ones=number_ones+1
        }
      }
      PC[j] = number_ones/rows
    }
    
    for (i in 1:rows){
      for (j in 1:columns){
        p = ( PR[i]+PC[j] )/2;
        r = runif(1) 
        if(r < p){  
          B[i,j] = 1;
        }
      }
      
    }
    
    # remove unconected species if present
    B = bipartite::empty(B)
    
    # skip if network has fewer than 2 rows or columns
    if(nrow(B)<2 || ncol(B)<2) {next}
    
    # compute nestedness
    null_nest[t] = compute_nestedness(B)
    
    # convert to igraph and calculate modularity
    graph = graph_from_incidence_matrix(B)
    modules = cluster_louvain(graph)
    null_mod[t] = modularity(modules)
    
    t=t+1
    
  }
  
  # unlist results 
  null_nest = unlist(null_nest)
  null_mod = unlist(null_mod)
  
  return(list(null_nest, null_mod))
  
}

compute_network_structure = function(network){
  
  # network incidence matrix
  M_inc = read.table(paste0("Data/",network,"/M_inc.csv"))
  M_inc = as.matrix(M_inc)
  
  # number of resources, consumers and interactions
  n_r = nrow(M_inc)
  n_c = ncol(M_inc)
  n_int = sum(M_inc)
  
  # connectance
  C = n_int / (n_r*n_c)
  
  # nestedness - observed
  nestedness_obs = compute_nestedness(M_inc)
  
  # modularity - observed
  graph = graph_from_incidence_matrix(M_inc)
  modules = cluster_louvain(graph)
  modularity_obs = modularity(modules)
  
  # run cell null model with current network and 100 replicates
  null_output = cell_model(M_inc,100)
  null_nest = null_output[[1]]
  null_mod = null_output[[2]]
  
  # compute mean of nestedness and modularity
  mean_nest = mean(null_nest, na.rm=TRUE)
  mean_mod = mean(null_mod, na.rm=TRUE)
  
  # compute standard devation of nestedness and modularity
  sd_nest = sd(null_nest, na.rm=TRUE)
  sd_mod = sd(null_mod, na.rm=TRUE)
  
  # compute zscore for nestedness and modularity
  z_score_nest = (nestedness_obs - mean_nest) / sd_nest
  z_score_mod = (modularity_obs - mean_mod) / sd_mod
  
  # store results
  df_out = data.frame(network=network,
                      n_res=n_r,
                      n_con=n_c,
                      con=C,
                      nest_obs=nestedness_obs,
                      nest_zscore=z_score_nest,
                      mod_obs=modularity_obs,
                      mod_zscore=z_score_mod)
  
  return(df_out)
}


# initialise dataframe for storing resutls
network_structure = data.frame(network=factor(),
                               n_con=integer(),
                               n_res=integer(),
                               con=double(),
                               nest_obs=double(),
                               nest_zscore=double(),
                               mod_obs=double(),
                               mod_zscore=double())

# compute network structure
for(network in networks){
  
  df_out = compute_network_structure(network)
  network_structure = rbind(network_structure, df_out)
  
}

# write out results
write.csv(network_structure, paste0("Results/networks_combined_M/metanetwork_structure.csv"), row.names=FALSE)
write.csv(network_structure, paste0("Results/networks_combined_A/metanetwork_structure.csv"), row.names=FALSE)


# PCA

network_structure = fread("Results/networks_combined_M/metanetwork_structure.csv")

network_structure = network_structure %>%
  mutate(n_sp=n_res+n_con)

pca_obs_nsp_con_nest_mod = PCA(network_structure %>% select(n_sp, con, nest_obs, mod_obs), graph=TRUE)
pca_obs_con_nest_mod = PCA(network_structure %>% select(con, nest_obs, mod_obs), graph=TRUE)

pca_zscore_nsp_con_nest_mod = PCA(network_structure %>% select(n_sp, con, nest_zscore, mod_zscore), graph=TRUE)
pca_zscore_con_nest_mod = PCA(network_structure %>% select(con, nest_zscore, mod_zscore), graph=TRUE)

network_structure_pca = network_structure %>%
  select(network, n_sp) %>%
  mutate(PC1_obs_nsp_con_nest_mod=as.vector(pca_obs_nsp_con_nest_mod$ind$coord[,1]),
         PC2_obs_nsp_con_nest_mod=as.vector(pca_obs_nsp_con_nest_mod$ind$coord[,2]),
         PC1_obs_con_nest_mod=as.vector(pca_obs_con_nest_mod$ind$coord[,1]))

# write out results
write.csv(network_structure_pca, paste0("Results/networks_combined_M/metanetwork_structure_pca.csv"), row.names=FALSE)
write.csv(network_structure_pca, paste0("Results/networks_combined_A/metanetwork_structure_pca.csv"), row.names=FALSE)
