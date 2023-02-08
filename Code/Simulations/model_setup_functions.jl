# FUNCTIONS USED TO SET UP SIMULATIONS

function get_network(networkName)

	# get network incidence matrix
	M_inc = readdlm(string("../../Data/",networkName,"/M_inc.csv"), ' ', Int)

	# number of resources
	n_r = size(M_inc, 1)

	# number of consumers
	n_c = size(M_inc, 2)

	# return
	return M_inc, n_r, n_c
end


function setup_habitat_destruction(dD, n)

	# number of patches
	n_patch = n * n

	# number of patches destroyed in a step
	dn = floor(Int, dD*n_patch)

	# vector of number of destroyed patches
	d = collect(0:dn:n_patch)

	# vector of fraction of destroyed patches
	D = d/n_patch

	# return
	return n_patch, dn, d, D
end


function setup_habitat_restoration(dD, n, DR)

	# number of patches
	n_patch = n * n

	# number of patches destroyed in a step
	dn = floor(Int, dD*n_patch)

	# vector of number of destroyed patches
	d = collect(DR*n_patch:-dn:0)

	# vector of fraction of destroyed patches
	D = d/n_patch

	# return
	return n_patch, dn, d, D
end


function setup_grids(n, n_r, n_c, d)

	# grid of patch states (1=pristine, 0=destroyed) for each D
	x_state = ones(Int, n, n, length(d))

	# grid of resources (1=present, 0=absent) for each species
	x_r = ones(Int, n, n, n_r)

	# grid of consumers (1=present, 0=absent) for each species
	x_c = ones(Int, n, n, n_c)

	# return
	return x_state, x_r, x_c
end


function transform_grid(n, n_patch, p_p)

	# destroyed patches / patches where species h is absent
	p_non_p = collect(1:n_patch)
	p_non_p = p_non_p[(!in).(p_non_p, Ref(p_p))]

	# combine present and absent patches
	p_all = vcat(hcat(p_p, ones(Int, length(p_p))), hcat(p_non_p, zeros(Int, length(p_non_p))))
	p_all = p_all[sortperm(p_all[:,1]), :]

	# reshape into nxn grid
	x_reshaped = reshape(p_all[:,2], (n,n))

	# return
	return x_reshaped
end


function setup_grids_restoration(n, n_patch, n_r, n_c, x_state, x_r, x_c, df_state, df_out, DR)

	# pristine patches
	p_p = df_state[df_state.D .=== DR,"patch_no"]

	# initial grid of patch states = grid from destruction at D=DR
	x_state[:,:,1] = transform_grid(n, n_patch, p_p)

	# initial grid of resources = grid from destruction at D=DR
	for h=1:n_r

		# patches where species h is present
		p_p = df_out[(df_out.guild .=== "resources") .& (df_out.species .=== h),"patch_no"]

		# reshape into nxn grid
		x_r[:,:,h] = transform_grid(n, n_patch, p_p)

	end

	# initial grid of consumers = grid from destruction at D=DR
	for h=1:n_c

		# patches where species h is present
		p_p = df_out[(df_out.guild .=== "consumers") .& (df_out.species .=== h),"patch_no"]

		# reshape into nxn grid
		x_c[:,:,h] = transform_grid(n, n_patch, p_p)

	end

	# return
	return x_state, x_r, x_c
end


function setup_moore_neighborhood(i, j, n)

	# Moore's neighborhood
	neigh = [i-1 j-1
             i+0 j-1
             i+1 j-1
             i-1 j+0
	         i+1 j+0
    	     i-1 j+1
             i+0 j+1
             i+1 j+1]

	# remove non-existent neighbours
	neigh = neigh[(neigh[:,1].>=1) .& (neigh[:,2].>=1) .& (neigh[:,1].<=n) .& (neigh[:,2].<=n),:]

	# return
	return neigh
end
