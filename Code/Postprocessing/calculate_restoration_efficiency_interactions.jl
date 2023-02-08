# CALCULATE INTERACTION RESTORATION EFFICIENCY

# this script calculates restoration efficiency for interactions

# Packages
@everywhere using DataFrames, CSV, DelimitedFiles

@everywhere function restoration_efficiency(networkName, networkType, dD, RS)

  # network incidence matrix
  M_inc = readdlm(string("../../Data/",networkName,"/M_inc.csv"), ' ', Int)

  # number of interactions
  n_int = sum(M_inc)

  # import abundance data
  abundanceD = CSV.read(string("../../Results/",networkName,"/",networkName,"_",networkType,"_D_abundance_interactions.csv"), DataFrame)
  abundanceR = CSV.read(string("../../Results/",networkName,"/",networkName,"_",networkType,"_R_",RS,"_abundance_interactions.csv"), DataFrame)

  # fractions of habitat destruction at which restoration starts
  DR_all = unique(abundanceR.DR)

  # create dataframe for storing results
  RE = DataFrame(DR=repeat(DR_all, outer=n_int),
                 interaction=repeat(1:n_int, inner=length(DR_all)),
                 RE=Vector{Union{Missing, Float16}}(missing, length(DR_all)*n_int))

  for DR in DR_all

    # habitat loss
    D = collect(DR:-dD:0)

    # abundance data for currrent DR
    abundanceD_DR = abundanceD[abundanceD.D.<=DR,:]
    abundanceR_DR = abundanceR[abundanceR.DR.===DR,:]

    for h=1:n_int

      # abundance - destruction
      abD = abundanceD_DR[abundanceD_DR.interaction.===h,"abundance"]

      # initialise area below destruction curve
      AD = 0

      # area below destruction curve
      for i=2:length(abD)
        AD = AD + 0.5*(abD[i-1]+abD[i])*dD
      end # i loop

      # abundance - restoration
      abR = abundanceR_DR[abundanceR_DR.interaction.===h,"abundance"]

      # initialise area below restoration curve
      AR = 0

      # area below restoration curve
      for i=2:length(abR)
        AR = AR + 0.5*(abR[i-1]+abR[i])*dD
      end # i loop

      # calculate restoration efficiency
      RE[(RE.DR.===DR) .& (RE.interaction.===h),"RE"] .= (AR - AD) / AD

    end # h loop - resources

  end # DR loop

  # Write out results
  CSV.write(string("../../Results/",networkName,"/",networkName,"_",networkType,"_",RS,"_restoration_efficiency_interactions.csv"), RE)

  return 0
end


# function for running networks in parallel
@everywhere function run_parallel(networkName, networkType, RS)

  dD = 0.05

  restoration_efficiency(networkName, networkType, dD, RS)

end


# Run simulations
#pmap(run_parallel, [networks], repeat([type], length(networks)), repeat([RS], length(networks)))

pmap(run_parallel, ["M_SD_025", "M_SD_027", "M_PL_036", "M_SD_008", "M_PL_059", "M_SD_014", "M_SD_005", "M_SD_002", "M_PL_033", "M_PL_037"], repeat(["M"], 10), repeat(["diversity"], 10))

#pmap(run_parallel, ["M_PL_025", "M_PL_046", "M_SD_010", "M_SD_012", "M_PL_039", "M_PL_006", "M_SD_007", "M_SD_016", "M_PL_051", "M_PL_010"], repeat(["M"], 10), repeat(["diversity"], 10))

#pmap(run_parallel, ["A_HP_015", "A_HP_035", "A_HP_028", "A_HP_005", "A_HP_021", "A_HP_007", "A_HP_032", "A_HP_012", "A_HP_008", "A_HP_047"], repeat(["A"], 10), repeat(["nonrandom"], 10))

#pmap(run_parallel, ["A_HP_002", "A_HP_033", "A_HP_042", "A_HP_046", "A_HP_025", "A_HP_050", "A_PH_007", "A_PH_004", "A_PH_005", "A_PH_006"], repeat(["A"], 10), repeat(["nonrandom"], 10))


#################
# TO RUN THIS SCRIPT
# REPLACE NETWORK NAMES, TYPES AND VERSION IN LINE 126
# DO THE FOLLOWING IN THE COMMAND LINE
# julia -p 10 calculate_restoration_efficiency_interactions.jl
