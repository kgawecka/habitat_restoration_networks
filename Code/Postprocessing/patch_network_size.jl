# PATCH NETWORK SIZE

# this script calculates species and interactions richness in patches

# Packages
@everywhere using DataFrames, DelimitedFiles, CSV


@everywhere function patch_interactions(n_int, intID, df_out)

  # create empty dataframe for storing results
  df_int = DataFrame(patch_no = Int[],
                     interaction = Int[])

  for i=1:n_int

    # interacting species
    resource_sp = intID[intID.interaction.===i,"resource"][1]
    consumer_sp = intID[intID.interaction.===i,"consumer"][1]

    # patches where resource or consumer present
    patches_resource = df_out[(df_out.guild.==="resources") .& (df_out.species.===resource_sp),"patch_no"]
    patches_consumer = df_out[(df_out.guild.==="consumers") .& (df_out.species.===consumer_sp),"patch_no"]

    # patches where both resource and consumer present
    patches_int = intersect(patches_resource, patches_consumer)

    # store patch numbers
    df_int_current = DataFrame(patch_no=patches_int, interaction=i)
    df_int = vcat(df_int, df_int_current)

  end

  return df_int
end


@everywhere function patch_size(networkName, networkType, n, dD, RS)

  # network incidence matrix
  M_inc = readdlm(string("../../Data/",networkName,"/M_inc.csv"), ' ', Int)

  # number of consumer and resource species, and interactions
  n_r = size(M_inc, 1)
  n_c = size(M_inc, 2)
  n_int = sum(M_inc)

  # interaction ID
  intID = CSV.read(string("../../Data/",networkName,"/interactionID.csv"), DataFrame)

  # fractions of habitat loss
  D_all = collect(0:dD:1)

  # fractions of habitat destruction at which restoration starts
  DR_all = collect(dD:dD:(1-dD))

  # create epmty dataframe for storing results
  df_size = DataFrame(DR = Float32[],
                      D = Float32[],
                      patch_no = Int[],
                      n_sp = Int[],
                      n_int = Int[])

  # DESTRUCTION
  for k=1:(length(D_all)-1)

    # Import simulation results
    df_out = CSV.read(string("../../Output/",networkName,"/",networkName,"_",networkType,"_D_out_D",floor(Int,D_all[k]*100),".csv"), DataFrame)

    if nrow(df_out) === 0
      break
    end

    # patch species richness
    df_out_g = groupby(df_out, :patch_no)
    no_sp = combine(df_out_g, nrow)
    rename!(no_sp, :nrow => :n_sp)

    # patch interaction presence
    df_int = patch_interactions(n_int, intID, df_out)

    # patch interaction richness
    df_int_g = groupby(df_int, :patch_no)
    no_int = combine(df_int_g, nrow)
    rename!(no_int, :nrow => :n_int)

    # combine species and interaction dataframes
    df_size_current = outerjoin(no_sp, no_int, on=:patch_no)
    df_size_current = hcat(DataFrame(DR=fill(D_all[k], nrow(df_size_current)),
                                     D=fill(D_all[k], nrow(df_size_current))), df_size_current)

    # Store in dataframe
    df_size = vcat(df_size, df_size_current)

  end # k loop

  # RESTORATION
  for d=1:length(DR_all)

    # check if results for DR exist
    if isfile(string("../../Output/",networkName,"/",networkName,"_",networkType,"_R_",RS,"_DR",floor(Int,DR_all[d]*100),"_out_D0.csv")) === false
      continue
    end

    # habitat loss
    D_all = collect(DR_all[d]:-dD:0)

    for k=2:length(D_all)

      # Import simulation results
      df_out = CSV.read(string("../../Output/",networkName,"/",networkName,"_",networkType,"_R_",RS,"_DR",floor(Int,DR_all[d]*100),"_out_D",floor(Int,D_all[k]*100),".csv"), DataFrame)

      if nrow(df_out) === 0
        break
      end

      # patch species richness
      df_out_g = groupby(df_out, :patch_no)
      no_sp = combine(df_out_g, nrow)
      rename!(no_sp, :nrow => :n_sp)

      # patch interaction presence
      df_int = patch_interactions(n_int, intID, df_out)

      # patch interaction richness
      df_int_g = groupby(df_int, :patch_no)
      no_int = combine(df_int_g, nrow)
      rename!(no_int, :nrow => :n_int)

      # combine species and interaction dataframes
      df_size_current = outerjoin(no_sp, no_int, on=:patch_no)
      df_size_current = hcat(DataFrame(DR=fill(DR_all[d], nrow(df_size_current)),
                                       D=fill(D_all[k], nrow(df_size_current))), df_size_current)

      # Store in dataframe
      df_size = vcat(df_size, df_size_current)

    end # k loop

  end # d loop

  # calculate fraction of species and interactions
  df_size.fract_sp = df_size.n_sp / (n_r+n_c)
  df_size.fract_int = df_size.n_int / n_int

  # Write out results
  CSV.write(string("../../Results/",networkName,"/",networkName,"_",networkType,"_R_",RS,"_patch_size.csv"), df_size)

  return 0
end


# function for running networks in parallel
@everywhere function run_parallel(networkName, networkType, RS)

  n = 100
  dD = 0.05

  patch_size(networkName, networkType, n, dD, RS)

end


# Run simulations
#pmap(run_parallel, [networks], repeat([type], length(networks)), repeat([RS], length(networks)))

pmap(run_parallel, ["M_SD_025", "M_SD_027", "M_PL_036", "M_SD_008", "M_PL_059", "M_SD_014", "M_SD_005", "M_SD_002", "M_PL_033", "M_PL_037"], repeat(["M"], 10), repeat(["diversity"], 10))

#pmap(run_parallel, ["M_PL_025", "M_PL_046", "M_SD_010", "M_SD_012", "M_PL_039", "M_PL_006", "M_SD_007", "M_SD_016", "M_PL_051", "M_PL_010"], repeat(["M"], 10), repeat(["diversity"], 10))

#pmap(run_parallel, ["A_HP_015", "A_HP_035", "A_HP_028", "A_HP_005", "A_HP_021", "A_HP_007", "A_HP_032", "A_HP_012", "A_HP_008", "A_HP_047"], repeat(["A"], 10), repeat(["nonrandom"], 10))

#pmap(run_parallel, ["A_HP_002", "A_HP_033", "A_HP_042", "A_HP_046", "A_HP_025", "A_HP_050", "A_PH_007", "A_PH_004", "A_PH_005", "A_PH_006"], repeat(["A"], 10), repeat(["nonrandom"], 10))


#################
# TO RUN THIS SCRIPT
# REPLACE NETWORK NAMES, TYPES AND VERSION IN LINE 184
# DO THE FOLLOWING IN THE COMMAND LINE
# julia -p 10 patch_network_size.jl