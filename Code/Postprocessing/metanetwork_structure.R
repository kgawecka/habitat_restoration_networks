# METANETWORK STRUCTURE

library(data.table)
library(dplyr)
library(igraph)
library(FactoMineR)

setwd("~/habitat_restoration_networks")


networks = c("C_0.1_R_1_S_20", "C_0.1_R_2_S_20", "C_0.1_R_3_S_20", "C_0.1_R_4_S_20", "C_0.1_R_5_S_20", 
             "C_0.2_R_1_S_20", "C_0.2_R_2_S_20", "C_0.2_R_3_S_20", "C_0.2_R_4_S_20", "C_0.2_R_5_S_20",
             "C_0.3_R_1_S_20", "C_0.3_R_2_S_20", "C_0.3_R_3_S_20", "C_0.3_R_4_S_20", "C_0.3_R_5_S_20", 
             "C_0.4_R_1_S_20", "C_0.4_R_2_S_20", "C_0.4_R_3_S_20", "C_0.4_R_4_S_20", "C_0.4_R_5_S_20",
             "C_0.1_R_1_S_40", "C_0.1_R_2_S_40", "C_0.1_R_3_S_40", "C_0.1_R_4_S_40", "C_0.1_R_5_S_40", 
             "C_0.2_R_1_S_40", "C_0.2_R_2_S_40", "C_0.2_R_3_S_40", "C_0.2_R_4_S_40", "C_0.2_R_5_S_40",
             "C_0.3_R_1_S_40", "C_0.3_R_2_S_40", "C_0.3_R_3_S_40", "C_0.3_R_4_S_40", "C_0.3_R_5_S_40", 
             "C_0.4_R_1_S_40", "C_0.4_R_2_S_40", "C_0.4_R_3_S_40", "C_0.4_R_4_S_40", "C_0.4_R_5_S_40",
             "C_0.1_R_1_S_80", "C_0.1_R_2_S_80", "C_0.1_R_3_S_80", "C_0.1_R_4_S_80", "C_0.1_R_5_S_80", 
             "C_0.2_R_1_S_80", "C_0.2_R_2_S_80", "C_0.2_R_3_S_80", "C_0.2_R_4_S_80", "C_0.2_R_5_S_80",
             "C_0.3_R_1_S_80", "C_0.3_R_2_S_80", "C_0.3_R_3_S_80", "C_0.3_R_4_S_80", "C_0.3_R_5_S_80", 
             "C_0.4_R_1_S_80", "C_0.4_R_2_S_80", "C_0.4_R_3_S_80", "C_0.4_R_4_S_80", "C_0.4_R_5_S_80")

networks = c("M_SD_005", "M_SD_008", "M_SD_010", "M_SD_012", "M_SD_025",
             "M_SD_002", "M_SD_007", "M_SD_014", "M_SD_016", "M_SD_027",
             "M_PL_006", "M_PL_010", "M_PL_036", "M_PL_037", "M_PL_059",
             "M_PL_025", "M_PL_033", "M_PL_039", "M_PL_046", "M_PL_051")

networks = c("A_HP_015", "A_HP_035", "A_HP_028", "A_HP_005", "A_HP_021", 
             "A_HP_007", "A_HP_032", "A_HP_012", "A_HP_008", "A_HP_047", 
             "A_HP_002", "A_HP_033", "A_HP_042", "A_HP_046", "A_HP_025",
             "A_HP_050", "A_PH_007", "A_PH_004", "A_PH_005", "A_PH_006")

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
write.csv(network_structure, paste0("Results/networks_combined/metanetwork_structure.csv"), row.names=FALSE)



# PCA

network_structure = fread("Results/networks_combined/metanetwork_structure.csv")

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
write.csv(network_structure_pca, paste0("Results/networks_combined/metanetwork_structure_pca.csv"), row.names=FALSE)


# ARTIFICIAL NETWORKS

# pca_obs_nsp_con_nest_mod

# PC1 = 71.45%; PC2 = 23.76%

#           Dim.1       Dim.2      Dim.3      Dim.4
#n_sp     -0.3381960  0.93983634 -0.0264060 0.04042071
#con       0.9704397 -0.02626039 -0.1697226 0.16956242
#nest_obs  0.9446082  0.15706580  0.2869680 0.02636245
#mod_obs  -0.9537441 -0.20442317  0.1208890 0.18430748

# pca_obs_con_nest_mod

# PC1 = 92.70%; PC2 = 4.45%
#           Dim.1       Dim.2      Dim.3
#con       0.9617516 -0.22345614 0.15843379
#nest_obs  0.9558978  0.28327253 0.07756327
#mod_obs  -0.9706861  0.05755753 0.23335710



# MUTUALISTIC NETWORKS

# pca_obs_nsp_con_nest_mod

# PC1 = 74.23%; PC2 = 18.71%

#            Dim.1        Dim.2       Dim.3       Dim.4
#n_sp     -0.6744115  0.729756262  0.06717529  0.09006877
#con       0.9703564 -0.005450768 -0.12094028  0.20917004
#nest_obs  0.9159063  0.124584952  0.38058389 -0.02738914
#mod_obs  -0.8567486 -0.447432369  0.21700671  0.13672644

# pca_obs_con_nest_mod

# PC1 = 87.23%; PC2 = 7.96%

#            Dim.1      Dim.2      Dim.3
#con       0.9514593 0.02296306  0.3069168
#nest_obs  0.9270152 0.33234186 -0.1737578
#mod_obs  -0.9231103 0.35741598  0.1418494



# ANTAGONISTIC NETWORKS

# pca_obs_nsp_con_nest_mod

# PC1 = 69.23%; PC2 = 21.45%

#           Dim.1       Dim.2      Dim.3       Dim.4
#n_sp     -0.5968142  0.79915919 0.03448746  0.06299189
#con       0.9866786 -0.07531291 0.01367332  0.14354897
#nest_obs  0.8606052  0.29457053 0.41007335 -0.06653338
#mod_obs  -0.8360077 -0.35615734 0.41365629  0.05596034

# pca_obs_con_nest_mod

# PC1 = 83.76%; PC2 = 11.46%

#             Dim.1      Dim.2      Dim.3
#con       0.9523098 0.04879799  0.3012054
#nest_obs  0.9044779 0.38226566 -0.1891896
#mod_obs  -0.8875237 0.44192807  0.1303888



# MUTUALISTIC & ANTAGONISTIC NETWORKS

# pca_obs_nsp_con_nest_mod

# PC1 = 70.84%; PC2 = 20.34%

#           Dim.1      Dim.2       Dim.3       Dim.4
#n_sp     -0.6183276  0.7791647  0.05962941  0.08377137
#con       0.9756552 -0.0387096 -0.06264138  0.20657820
#nest_obs  0.8838711  0.1966060  0.41942590 -0.06480633
#mod_obs  -0.8473278 -0.4080734  0.32187251  0.10913190

# pca_obs_con_nest_mod

# PC1 = 85.03%; PC2 = 9.96%

#           Dim.1      Dim.2      Dim.3
#con       0.9502441 0.01061917  0.3113252
#nest_obs  0.9091439 0.38029345 -0.1698067
#mod_obs  -0.9063761 0.39258788  0.1560679
