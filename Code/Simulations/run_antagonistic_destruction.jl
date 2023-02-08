# load packages
using Distributed, DataFrames, Random, CSV, DelimitedFiles, StatsBase
@everywhere using Distributed, DataFrames, Random, CSV, DelimitedFiles, StatsBase

# load functions
@everywhere include("model_setup_functions.jl")
@everywhere include("master_functions.jl")
@everywhere include("extinction_colonisation_A_functions.jl")
@everywhere include("output_A_D_functions.jl")
@everywhere include("habitat_destruction_functions.jl")


# define model function
@everywhere function dynamicsD(networkName, n, dD, tmin, tmax, tol, er, ec, cr, cc)

  # get network incidence matrix
  M_inc, n_r, n_c = get_network(networkName)

  # setup variables for habitat destruction
  n_patch, dn, d, D = setup_habitat_destruction(dD, n)

  # setup grids
  x_state, x_r, x_c = setup_grids(n, n_r, n_c, d)

  # loop through habitat destruction fractions
  @inbounds for k=1:(length(d)-1)

    # run colonisations/extinctions until steady stated
    x_r, x_c = iterate_model_until_steady_state(x_state, x_r, x_c, n_r, n_c, M_inc, er, ec, cr, cc, n, n_patch, tmax, tmin, tol, k)

    # store and output results
    store_species_results = store_and_write_species_results(networkName, n_r, n_c, x_r, x_c, n_patch, D, k)

    # habitat destruction
    x_state, x_r, x_c = habitat_destruction(n, n_patch, dn, x_state, x_r, x_c, k)

  end # k loop

  # store and output results
  store_patch_results = store_and_write_patch_results(networkName, x_state, n_patch, D)

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
  er = 0.10
  ec = 0.40
  cr = 0.10
  cc = 0.10

  dynamicsD(networkName, n, dD, tmin, tmax, tol, er, ec, cr, cc)

end

# Run simulations
#pmap(run_parallel, [networkName])

#pmap(run_parallel, ["C_0.1_R_1_S_40", "C_0.1_R_2_S_40", "C_0.1_R_3_S_40", "C_0.1_R_4_S_40", "C_0.1_R_5_S_40", "C_0.2_R_1_S_40", "C_0.2_R_2_S_40", "C_0.2_R_3_S_40", "C_0.2_R_4_S_40", "C_0.2_R_5_S_40"])

#pmap(run_parallel, ["C_0.3_R_1_S_40", "C_0.3_R_2_S_40", "C_0.3_R_3_S_40", "C_0.3_R_4_S_40", "C_0.3_R_5_S_40", "C_0.4_R_1_S_40", "C_0.4_R_2_S_40", "C_0.4_R_3_S_40", "C_0.4_R_4_S_40", "C_0.4_R_5_S_40"])

#pmap(run_parallel, ["A_HP_015", "A_HP_035", "A_HP_028", "A_HP_005", "A_HP_021", "A_HP_007", "A_HP_032", "A_HP_012", "A_HP_008", "A_HP_047"])

#pmap(run_parallel, ["A_HP_002", "A_HP_033", "A_HP_042", "A_HP_046", "A_HP_025", "A_HP_050", "A_PH_007", "A_PH_004", "A_PH_005", "A_PH_006"])


#################
# TO RUN THIS SCRIPT
# DO THE FOLLOWING IN THE COMMAND LINE, REPLACE "5" WITH NUMBER OF CORES (NETWORKS)
# julia -p 10 run_antagonistic_destruction.jl
