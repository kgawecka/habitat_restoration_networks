# CALCULATE INTERACTION ABUNDANCE

# this script calculates interaction abundance in habitat destruction or restoration simulations

# Packages
@everywhere using DataFrames, DelimitedFiles, CSV


@everywhere function abundance_destruction(networkName, networkType, n, dD)

  # network incidence matrix
  M_inc = readdlm(string("../../Data/",networkName,"/M_inc.csv"), ' ', Int)

  # number of consumer and resource species, and interactions
  n_r = size(M_inc, 1)
  n_c = size(M_inc, 2)
  n_int = sum(M_inc)

  # interaction ID
  intID = CSV.read(string("../../Data/",networkName,"/interactionID.csv"), DataFrame)

  # habitat loss
  D = collect(0:dD:1)

  # create empty dataframe
  ab_df = DataFrame(D = Float32[],
                    interaction = Int[],
                    abundance = Float32[])

  for k=1:(length(D)-1)

    # Import simulation results
    df_out = CSV.read(string("../../Output/",networkName,"/",networkName,"_",networkType,"_D_out_D",floor(Int,D[k]*100),".csv"), DataFrame)

    if nrow(df_out) === 0
      break
    end

    # create empty vector for storing interaction abunndances
    int_ab = zeros(Float16, n_int)

    for i=1:n_int

      # interacting species
      resource_sp = intID[intID.interaction.===i,"resource"][1]
      consumer_sp = intID[intID.interaction.===i,"consumer"][1]

      # patches where resource or consumer present
      patches_resource = df_out[(df_out.guild.==="resources") .& (df_out.species.===resource_sp),"patch_no"]
      patches_consumer = df_out[(df_out.guild.==="consumers") .& (df_out.species.===consumer_sp),"patch_no"]

      # patches where both resource and consumer present
      patches_int = intersect(patches_resource, patches_consumer)

      # store number of patches
      int_ab[i] = length(patches_int) / (n*n)

    end

    # Store in dataframe
    ab = DataFrame(D=fill(D[k], n_int), interaction=1:n_int, abundance=int_ab)
    ab_df = vcat(ab_df, ab)

  end

  # Create dataframe with all D and species
  abundance = DataFrame(D = repeat(D, inner=n_int),
                        interaction = repeat(1:n_int, outer=length(D)))
  abundance = leftjoin(abundance, ab_df, on=[:D, :interaction])
  abundance.abundance = coalesce.(abundance.abundance, 0.0)

  # Write out results
  CSV.write(string("../../Results/",networkName,"/",networkName,"_",networkType,"_D_abundance_interactions.csv"), abundance)

  return 0
end


@everywhere function abundance_restoration(networkName, networkType, n, dD, RS)

  # network incidence matrix
  M_inc = readdlm(string("../../Data/",networkName,"/M_inc.csv"), ' ', Int)

  # number of consumer and resource species, and interactions
  n_r = size(M_inc, 1)
  n_c = size(M_inc, 2)
  n_int = sum(M_inc)

  # interaction ID
  intID = CSV.read(string("../../Data/",networkName,"/interactionID.csv"), DataFrame)

  # fractions of habitat destruction at which restoration starts
  DR_all = collect(dD:dD:(1-dD))

  # create empty dataframe
  ab_df = DataFrame(DR = Float32[],
                    D = Float32[],
                    interaction = Int[],
                    abundance = Float32[])

  for d=1:length(DR_all)

    # check if results for DR exist
    if isfile(string("../../Output/",networkName,"/",networkName,"_",networkType,"_R_",RS,"_DR",floor(Int,DR_all[d]*100),"_out_D0.csv")) === false
      continue
    end

    # habitat loss
    D = collect(DR_all[d]:-dD:0)

    for k=1:length(D)

      # Import simulation results
      if k === 1
        df_out = CSV.read(string("../../Output/",networkName,"/",networkName,"_",networkType,"_D_out_D",floor(Int,DR_all[d]*100),".csv"), DataFrame)
      else
        df_out = CSV.read(string("../../Output/",networkName,"/",networkName,"_",networkType,"_R_",RS,"_DR",floor(Int,DR_all[d]*100),"_out_D",floor(Int,D[k]*100),".csv"), DataFrame)
      end

      if nrow(df_out) === 0
        break
      end

      # create empty vector for storing interaction abunndances
      int_ab = zeros(Float16, n_int)

      for i=1:n_int

        # interacting species
        resource_sp = intID[intID.interaction.===i,"resource"][1]
        consumer_sp = intID[intID.interaction.===i,"consumer"][1]

        # patches where resource or consumer present
        patches_resource = df_out[(df_out.guild.==="resources") .& (df_out.species.===resource_sp),"patch_no"]
        patches_consumer = df_out[(df_out.guild.==="consumers") .& (df_out.species.===consumer_sp),"patch_no"]

        # patches where both resource and consumer present
        patches_int = intersect(patches_resource, patches_consumer)

        # store number of patches
        int_ab[i] = length(patches_int) / (n*n)

      end # i loop

      # Store in dataframe
      ab = DataFrame(DR=fill(DR_all[d], n_int),
                     D=fill(D[k], n_int),
                     interaction=1:n_int,
                     abundance=int_ab)
      ab_df = vcat(ab_df, ab)

    end # k loop
  end # d loop

  # Create dataframe with all D and species
  DRs = unique(ab_df.DR)
  D = collect(0:dD:1)
  abundance = DataFrame(DR = repeat(DRs, inner=(length(D)*n_int)),
                        D = repeat(D, outer=length(DRs), inner=n_int),
                        interaction = repeat(1:n_int, outer=(length(DRs)*length(D))))
  abundance = abundance[abundance.D.<=abundance.DR,:]
  abundance = leftjoin(abundance, ab_df, on=[:DR, :D, :interaction])
  abundance.abundance = coalesce.(abundance.abundance, 0.0)

  # Write out results
  CSV.write(string("../../Results/",networkName,"/",networkName,"_",networkType,"_R_",RS,"_abundance_interactions.csv"), abundance)

  return 0
end


# function for running networks in parallel
@everywhere function run_parallel(networkName, networkType, habitatChange)

  n = 100
  dD = 0.05

  if habitatChange === "destruction"

    abundance_destruction(networkName, networkType, n, dD)

  elseif habitatChange === "restoration"

    RS = "diversity"

    abundance_restoration(networkName, networkType, n, dD, RS)

  end

end


# Run simulations
#pmap(run_parallel, [networks], repeat([type], length(networks)), repeat([habitatChange], length(networks)))

pmap(run_parallel, ["M_SD_025", "M_SD_027", "M_PL_036", "M_SD_008", "M_PL_059", "M_SD_014", "M_SD_005", "M_SD_002", "M_PL_033", "M_PL_037"], repeat(["M"], 10), repeat(["restoration"], 10))

#pmap(run_parallel, ["M_PL_025", "M_PL_046", "M_SD_010", "M_SD_012", "M_PL_039", "M_PL_006", "M_SD_007", "M_SD_016", "M_PL_051", "M_PL_010"], repeat(["M"], 10), repeat(["restoration"], 10))

#pmap(run_parallel, ["A_HP_015", "A_HP_035", "A_HP_028", "A_HP_005", "A_HP_021", "A_HP_007", "A_HP_032", "A_HP_012", "A_HP_008", "A_HP_047"], repeat(["A"], 10), repeat(["restoration"], 10))

#pmap(run_parallel, ["A_HP_002", "A_HP_033", "A_HP_042", "A_HP_046", "A_HP_025", "A_HP_050", "A_PH_007", "A_PH_004", "A_PH_005", "A_PH_006"], repeat(["A"], 10), repeat(["restoration"], 10))


#################
# TO RUN THIS SCRIPT
# REPLACE NETWORK NAMES, TYPES AND VERSION IN LINE 184
# DO THE FOLLOWING IN THE COMMAND LINE
# julia -p 10 calculate_abundance_interactions.jl
