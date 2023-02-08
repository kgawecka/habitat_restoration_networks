# load packages
using Distributed, DataFrames, Random, CSV, DelimitedFiles, StatsBase
@everywhere using Distributed, DataFrames, Random, CSV, DelimitedFiles, StatsBase

# load functions
@everywhere include("model_setup_functions.jl")
@everywhere include("master_functions.jl")
@everywhere include("extinction_colonisation_M_functions.jl")
@everywhere include("output_M_R_functions.jl")
@everywhere include("habitat_restoration_functions.jl")


# define model function
@everywhere function dynamicsD(networkName, n, dD, tmin, tmax, tol, er, ec, cr, cc, patch_neigh, DR_all, RS)

  # get network incidence matrix
  M_inc, n_r, n_c = get_network(networkName)
  
  if RS === "generalist"
    # find most generalist species
    target_guild = ifelse(sum(M_inc[1,:]) >= sum(M_inc[:,1]), "resource", "consumer")
    target_species = 1
  elseif RS === "specialist"
    # find most specialist species
    target_guild = ifelse(sum(M_inc[n_r,:]) <= sum(M_inc[:,n_c]), "resource", "consumer")
    target_species = ifelse(target_guild === "resource", n_r, n_c)
  end

  # patch state from destruction simulations
  df_state = CSV.read(string("../../Output/",networkName,"/",networkName,"_M_D_patch_state.csv"), DataFrame)

  # loop through starting habitat restoration fractions
  @inbounds for r=1:length(DR_all)

    # current DR
    DR = DR_all[r]

    # initial species grid for restoration
  	df_out = CSV.read(string("../../Output/",networkName,"/",networkName,"_M_D_out_D",floor(Int,DR*100),".csv"), DataFrame)

    # if no species present in landscape, move to next DR
    if nrow(df_out) === 0
      continue
    end

    # setup variables for habitat restoration
    n_patch, dn, d, D = setup_habitat_restoration(dD, n, DR)

    # setup grids
    x_state, x_r, x_c = setup_grids(n, n_r, n_c, d)

    # initial grids for restoration
    x_state, x_r, x_c = setup_grids_restoration(n, n_patch, n_r, n_c, x_state, x_r, x_c, df_state, df_out, DR)

    # loop through habitat destruction fractions
    @inbounds for k=2:(length(d))

      # habitat restoration
      if RS === "random"
        x_state = habitat_restoration_random(n, n_patch, dn, x_state, k)
      elseif RS === "nonrandom"
        x_state = habitat_restoration_nonrandom(n, n_patch, dn, patch_neigh, x_state, x_r, x_c, k)
      elseif RS === "diversity"
        x_state = habitat_restoration_target_diversity(n, n_patch, dn, patch_neigh, x_state, x_r, x_c, k)
      else
        x_state = habitat_restoration_target_species(n, n_patch, dn, patch_neigh, x_state, x_r, x_c, target_guild, target_species, k)
      end

      # run colonisations/extinctions until steady stated
      x_r, x_c = iterate_model_until_steady_state(x_state, x_r, x_c, n_r, n_c, M_inc, er, ec, cr, cc, n, n_patch, tmax, tmin, tol, k)

      # store and output results
      store_species_results = store_and_write_species_results(networkName, n_r, n_c, x_r, x_c, n_patch, D, k, DR, RS)

    end # k loop

    # store and output results
    store_patch_results = store_and_write_patch_results(networkName, x_state, n_patch, D, DR, RS)

  end # r loop

# return
return 0
end


# define function for running in parallel
@everywhere function run_parallel(networkName)

  # simulation parameters
  n = 100            # grid size nxn
  dD = 0.05          # fraction of patches destoryed at a time
  tmin = 10          # minimum number of timesteps
  tmax = 1000        # maximum number of timesteps
  tol = 0.001        # convergence tolerance

  # extinction and colonisation probabilities
  er = 0.20
  ec = 0.20
  cr = 0.10
  cc = 0.10

  # fractions of habitat destruction at which restoration starts
  #DR_all = [0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95]
  DR_all = collect(dD:dD:(1-dD))

  # restoration strategy: random, nonrandom, diversity, generalist, specialist
  RS = "diversity"

  # patch neighbourhood
  patch_neigh = CSV.read("../../Data/patch_neighbourhood.csv", DataFrame)

  dynamicsD(networkName, n, dD, tmin, tmax, tol, er, ec, cr, cc, patch_neigh, DR_all, RS)

end

# Run simulations
#pmap(run_parallel, [networkName])

#pmap(run_parallel, ["C_0.1_R_1_S_80", "C_0.1_R_2_S_80", "C_0.1_R_3_S_80", "C_0.1_R_4_S_80", "C_0.1_R_5_S_80", "C_0.2_R_1_S_80", "C_0.2_R_2_S_80", "C_0.2_R_3_S_80", "C_0.2_R_4_S_80", "C_0.2_R_5_S_80"])
#pmap(run_parallel, ["C_0.3_R_1_S_80", "C_0.3_R_2_S_80", "C_0.3_R_3_S_80", "C_0.3_R_4_S_80", "C_0.3_R_5_S_80", "C_0.4_R_1_S_80", "C_0.4_R_2_S_80", "C_0.4_R_3_S_80", "C_0.4_R_4_S_80", "C_0.4_R_5_S_80"])

pmap(run_parallel, ["M_SD_025", "M_SD_027", "M_PL_036", "M_SD_008", "M_PL_059", "M_SD_014", "M_SD_005", "M_SD_002", "M_PL_033", "M_PL_037"])

#pmap(run_parallel, ["M_PL_025", "M_PL_046", "M_SD_010", "M_SD_012", "M_PL_039", "M_PL_006", "M_SD_007", "M_SD_016", "M_PL_051", "M_PL_010"])


#################
# TO RUN THIS SCRIPT
# DO THE FOLLOWING IN THE COMMAND LINE, REPLACE "5" WITH NUMBER OF CORES (NETWORKS)
# julia -p 10 run_mutualistic_restoration.jl
