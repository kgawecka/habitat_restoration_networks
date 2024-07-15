# Habitat restoration and the recovery of metacommunities

This repository contains the code used to perform simulations and postprocess the results for the manuscript "Habitat restoration and the recovery of metacommunities" [Gawecka & Bascompte, 2023](https://doi.org/10.1111/1365-2664.14445).

## `Code`
All code was created in either Julia version 1.4.2 or R version 3.6.2.

The workflow is structured as follows:

1. Run habitat destruction simulations with spatially explicit metacommunity models
    - `Simulations/run_antagonistic_destruction.jl` - antagonistic networks
    - `Simulations/run_mutualistic_destruction.jl` - mutualistic networks

2. Run habitat restoration simulations with spatially explicit metacommunity models
    - `Simulations/run_antagonistic_restoration.jl` - antagonistic networks
    - `Simulations/run_mutualistic_restoration.jl` - mutualistic networks

3. Postprocess results
    - `Postprocessing/calculate_abundance_interactions.jl` - computes regional interaction abundances
    - `Postprocessing/patch_network_size.jl` - computes size of local networks
    - `Postprocessing/patch_dissimilarity.R` - computes beta-diversity between pairs of patches and with the metanetwork
    - `Postprocessing/calculate_restoration_efficiency_interactions.jl` - computes restoration efficiency of regional interaction abundances
    - `Postprocessing/calculate_restoration_efficiency_networks.jl` - computes restoration efficiency of local networks (size and beta-diversity)
    - `Postprocessing/metanetwork_structure.R` - PCA analysis of structure of metanetworks

4. Combine results of all networks
    - `Postprocessing/combine_networks_M.R` - antagonistic networks
    - `Postprocessing/combine_networks_A.R` - mutualistic networks

## `Data`
This directory contains the incidence matrices and interaction ID for all networks used in the simulations for the manuscript. It is divided into directories named after each network.

## `Output` and `Results`
When running simulations and postprocessing, `Output` and `Results` directories must be created. 
Additionally, these directories must contin subdirectories named after the network being simulated (as in `Data` directory). 
