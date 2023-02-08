# CALCULATE NETWORK RESTORATION EFFICIENCY

# this script calculates restoration efficiency for network size and beta diversity

# Packages
@everywhere using DataFrames, CSV, DelimitedFiles, Statistics

@everywhere function calculate_R(dataD, dataR, dD)

  # initialise area below destruction curve
  AD = 0

  # area below destruction curve
  for i=2:length(dataD)
    AD = AD + 0.5*(dataD[i-1]+dataD[i])*dD
  end # i loop

  # initialise area below restoration curve
  AR = 0

  # area below restoration curve
  for i=2:length(dataR)
    AR = AR + 0.5*(dataR[i-1]+dataR[i])*dD
  end # i loop

  # calculate restoration efficiency
  R = (AR - AD) / AD

  return R
end

@everywhere function restoration_efficiency(networkName, networkType, dD, RS)

  # import patch network size and beta diversity data
  patch_size = CSV.read(string("../../Results/",networkName,"/",networkName,"_",networkType,"_R_",RS,"_patch_size.csv"), DataFrame)
  patch_betadiv = CSV.read(string("../../Results/",networkName,"/",networkName,"_",networkType,"_R_",RS,"_patch_dissimilarity.csv"), DataFrame)

  # calculate mean values across all patches
  patch_size_sp = combine(groupby(patch_size, [:DR, :D]), :fract_sp => mean∘skipmissing)
  patch_size_int = combine(groupby(patch_size, [:DR, :D]), :fract_int => mean∘skipmissing)
  betadiv_patch_sp = combine(groupby(patch_betadiv[patch_betadiv.patch_j.>0,:], [:DR, :D]), :beta_S => mean∘skipmissing)
  betadiv_patch_int = combine(groupby(patch_betadiv[patch_betadiv.patch_j.>0,:], [:DR, :D]), :beta_WN => mean∘skipmissing)
  betadiv_meta_sp = combine(groupby(patch_betadiv[patch_betadiv.patch_j.===0,:], [:DR, :D]), :beta_S => mean∘skipmissing)
  betadiv_meta_int = combine(groupby(patch_betadiv[patch_betadiv.patch_j.===0,:], [:DR, :D]), :beta_WN => mean∘skipmissing)
  replace!(patch_size_sp.fract_sp_mean_skipmissing, NaN => 0)
  replace!(patch_size_int.fract_int_mean_skipmissing, NaN => 0)
  replace!(betadiv_patch_sp.beta_S_mean_skipmissing, NaN => 0)
  replace!(betadiv_patch_int.beta_WN_mean_skipmissing, NaN => 0)
  replace!(betadiv_meta_sp.beta_S_mean_skipmissing, NaN => 0)
  replace!(betadiv_meta_int.beta_WN_mean_skipmissing, NaN => 0)

  # fractions of habitat destruction at which restoration starts
  DR_all = unique(vcat(patch_size_sp.DR, patch_size_int.DR, betadiv_patch_sp.DR, betadiv_patch_int.DR, betadiv_meta_sp.DR, betadiv_meta_int.DR))

  # create dataframe for storing results
  RE = DataFrame(DR=DR_all,
                 RE_size_sp=Vector{Union{Missing, Float16}}(missing, length(DR_all)),
                 RE_size_int=Vector{Union{Missing, Float16}}(missing, length(DR_all)),
                 RE_beta_patch_sp=Vector{Union{Missing, Float16}}(missing, length(DR_all)),
                 RE_beta_patch_int=Vector{Union{Missing, Float16}}(missing, length(DR_all)),
                 RE_beta_meta_sp=Vector{Union{Missing, Float16}}(missing, length(DR_all)),
                 RE_beta_meta_int=Vector{Union{Missing, Float16}}(missing, length(DR_all)))

  for d=2:length(DR_all)

    DR_current = DR_all[d]

    # habitat loss
    D_all = collect(DR_current:-dD:0)
    

    # NETWORK SIZE - SPECIES

    # network size - destruction
    dataD = patch_size_sp[(patch_size_sp.DR.==patch_size_sp.D) .& (patch_size_sp.D.<=DR_current),:]
    dataD = sort!(dataD, :D).fract_sp_mean_skipmissing

    # network size - restoration
    dataR = patch_size_sp[patch_size_sp.DR.===DR_current,:]
    dataR = sort!(dataR, :D).fract_sp_mean_skipmissing
    if length(dataR) < length(dataD)
      dataR = vcat(zeros(length(dataD)-length(dataR)), dataR)
    end

    # calculate restoration efficiency
    RE[RE.DR.===DR_current,"RE_size_sp"] .= calculate_R(dataD, dataR, dD)
    
    
    # NETWORK SIZE - INTERACTIONS

    # network size - destruction
    dataD = patch_size_int[(patch_size_int.DR.==patch_size_int.D) .& (patch_size_int.D.<=DR_current),:]
    dataD = sort!(dataD, :D).fract_int_mean_skipmissing

    # network size - restoration
    dataR = patch_size_int[patch_size_int.DR.===DR_current,:]
    dataR = sort!(dataR, :D).fract_int_mean_skipmissing
    if length(dataR) < length(dataD)
      dataR = vcat(zeros(length(dataD)-length(dataR)), dataR)
    end

    # calculate restoration efficiency
    RE[RE.DR.===DR_current,"RE_size_int"] .= calculate_R(dataD, dataR, dD)
    

    # BETA DIVERSITY - PATCHES - SPECIES

    # betadiv patches - destruction
    dataD = betadiv_patch_sp[(betadiv_patch_sp.DR.==betadiv_patch_sp.D) .& (betadiv_patch_sp.D.<=DR_current),:]
    dataD = sort!(dataD, :D).beta_S_mean_skipmissing
    dataD = 1 .- dataD

    # betadiv patches - restoration
    dataR = betadiv_patch_sp[betadiv_patch_sp.DR.===DR_current,:]
    dataR = sort!(dataR, :D).beta_S_mean_skipmissing
    if length(dataR) < length(dataD)
      dataR = vcat(ones(length(dataD)-length(dataR)), dataR)
    end
    dataR = 1 .- dataR

    # calculate restoration efficiency
    RE[RE.DR.===DR_current,"RE_beta_patch_sp"] .= calculate_R(dataD, dataR, dD)
    
    
    # BETA DIVERSITY - PATCHES - INTERACTIONS

    # betadiv patches - destruction
    dataD = betadiv_patch_int[(betadiv_patch_int.DR.==betadiv_patch_int.D) .& (betadiv_patch_int.D.<=DR_current),:]
    dataD = sort!(dataD, :D).beta_WN_mean_skipmissing
    dataD = 1 .- dataD

    # betadiv patches - restoration
    dataR = betadiv_patch_int[betadiv_patch_int.DR.===DR_current,:]
    dataR = sort!(dataR, :D).beta_WN_mean_skipmissing
    if length(dataR) < length(dataD)
      dataR = vcat(ones(length(dataD)-length(dataR)), dataR)
    end
    dataR = 1 .- dataR

    # calculate restoration efficiency
    RE[RE.DR.===DR_current,"RE_beta_patch_int"] .= calculate_R(dataD, dataR, dD)
    

    # BETA DIVERSITY - METAWEB - SPECIES

    # betadiv metaweb - destruction
    dataD = betadiv_meta_sp[(betadiv_meta_sp.DR.==betadiv_meta_sp.D) .& (betadiv_meta_sp.D.<=DR_current),:]
    dataD = sort!(dataD, :D).beta_S_mean_skipmissing
    dataD = 1 .- dataD

    # betadiv metaweb - restoration
    dataR = betadiv_meta_sp[betadiv_meta_sp.DR.===DR_current,:]
    dataR = sort!(dataR, :D).beta_S_mean_skipmissing
    if length(dataR) < length(dataD)
      dataR = vcat(ones(length(dataD)-length(dataR)), dataR)
    end
    dataR = 1 .- dataR

    # calculate restoration efficiency
    RE[RE.DR.===DR_current,"RE_beta_meta_sp"] .= calculate_R(dataD, dataR, dD)
    
    
    # BETA DIVERSITY - METAWEB - INTERACTIONS

    # betadiv metaweb - destruction
    dataD = betadiv_meta_int[(betadiv_meta_int.DR.==betadiv_meta_int.D) .& (betadiv_meta_int.D.<=DR_current),:]
    dataD = sort!(dataD, :D).beta_WN_mean_skipmissing
    dataD = 1 .- dataD

    # betadiv metaweb - restoration
    dataR = betadiv_meta_int[betadiv_meta_int.DR.===DR_current,:]
    dataR = sort!(dataR, :D).beta_WN_mean_skipmissing
    if length(dataR) < length(dataD)
      dataR = vcat(ones(length(dataD)-length(dataR)), dataR)
    end
    dataR = 1 .- dataR

    # calculate restoration efficiency
    RE[RE.DR.===DR_current,"RE_beta_meta_int"] .= calculate_R(dataD, dataR, dD)


  end # DR loop

  # Write out results
  CSV.write(string("../../Results/",networkName,"/",networkName,"_",networkType,"_",RS,"_restoration_efficiency_networks.csv"), RE)

  return 0
end


# function for running networks in parallel
@everywhere function run_parallel(networkName, networkType, RS)

  dD = 0.05

  restoration_efficiency(networkName, networkType, dD, RS)

end


# Run simulations
#pmap(run_parallel, [networks], repeat([type], length(networks)), repeat([RS], length(networks)))

#pmap(run_parallel, ["M_SD_025", "M_SD_027", "M_PL_036", "M_SD_008", "M_PL_059", "M_SD_014", "M_SD_005", "M_SD_002", "M_PL_033", "M_PL_037"], repeat(["M"], 10), repeat(["diversity"], 10))

#pmap(run_parallel, ["M_PL_025", "M_PL_046", "M_SD_010", "M_SD_012", "M_PL_039", "M_PL_006", "M_SD_007", "M_SD_016", "M_PL_051", "M_PL_010"], repeat(["M"], 10), repeat(["diversity"], 10))

#pmap(run_parallel, ["A_HP_015", "A_HP_035", "A_HP_028", "A_HP_005", "A_HP_021", "A_HP_007", "A_HP_032", "A_HP_012", "A_HP_008", "A_HP_047"], repeat(["A"], 10), repeat(["random"], 10))

#pmap(run_parallel, ["A_HP_002", "A_HP_033", "A_HP_042", "A_HP_046", "A_HP_025", "A_HP_050", "A_PH_007", "A_PH_004", "A_PH_005", "A_PH_006"], repeat(["A"], 10), repeat(["random"], 10))


#################
# TO RUN THIS SCRIPT
# REPLACE NETWORK NAMES, TYPES AND VERSION IN LINE 126
# DO THE FOLLOWING IN THE COMMAND LINE
# julia -p 10 calculate_restoration_efficiency_networks.jl
