# FUNCTIONS USED TO STORE AND WITE OUT RESULTS

function build_output_data(n_r, n_c, x_r, x_c, n_patch)

	# create dataframe
	df_out = vcat(DataFrame(patch_no = repeat(1:n_patch, outer=n_r),
												  guild = "resources",
									 			  species = repeat(1:n_r, inner=n_patch),
                   			  p = reshape(x_r, (n_patch*n_r))),
							  DataFrame(patch_no = repeat(1:n_patch, outer=n_c),
												  guild = "consumers",
												  species = repeat(1:n_c, inner=n_patch),
											    p = reshape(x_c, (n_patch*n_c))))

	# store species presence only
	df_out = df_out[df_out.p.===1,["patch_no", "guild", "species"]]

	# return
	return df_out
end


function store_and_write_species_results(networkName, n_r, n_c, x_r, x_c, n_patch, D, k, DR, RS)

	# store results
	df_out = build_output_data(n_r, n_c, x_r, x_c, n_patch)

	# output results
	CSV.write(string("../../Output/",networkName,"/",networkName,"_M_R_",RS,"_DR",floor(Int,DR*100),"_out_D",floor(Int,D[k]*100),".csv"), df_out)

	# return
	return 1
end


function store_and_write_patch_results(networkName, x_state, n_patch, D, DR, RS)

	# store results
	df_state = DataFrame(D = repeat(D, inner=n_patch),
						 patch_no = repeat(1:n_patch, outer=length(D)),
						 state = vec(x_state))

	# store pristine patches only
	df_state = df_state[df_state.state.===1,["D", "patch_no"]]

	# output results
	CSV.write(string("../../Output/",networkName,"/",networkName,"_M_R_",RS,"_DR",floor(Int,DR*100),"_patch_state.csv"), df_state)

	# return
	return 1
end
