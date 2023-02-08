# FUNCTIONS USED TO ITERATE THROUGH GRID AND TIMESTEPS

function update_patch_sequentially(M_inc, n, n_r, n_c, x_state, x_r, x_c, x_r_old, x_c_old, er, ec, cr, cc, k)

	# update each patch sequentially
	for j=1:n
		for i=1:n

			# check state of the current patch
			if x_state[i,j,k] === 0
				# if patch is destroyed, move to next patch
				continue
			end

			# setup moore neighborhood
			neigh = setup_moore_neighborhood(i, j, n)

			# resource extinctions and colonisations
			x_r = resource_extinctions_and_colonisations(M_inc, er, cr, n_r, x_r_old, x_r, x_c_old, neigh, i, j)

			# consumer extinctions and colonisations
			x_c = consumer_extinctions_and_colonisations(M_inc, ec, cc, n_c, x_c_old, x_c, x_r_old, neigh, i, j)

		end # i loop
	end # j loop

	# return
	return x_r, x_c
end


function check_model_convergence(n_patch, n_r, n_c, x_r_old, x_c_old, x_r, x_c)

	  p_r_old = sum(reshape(x_r_old, (n_patch, n_r)), dims=1) ./ n_patch
    p_r_new = sum(reshape(x_r, (n_patch, n_r)), dims=1) ./ n_patch
    p_c_old = sum(reshape(x_c_old, (n_patch, n_c)), dims=1) ./ n_patch
	  p_c_new = sum(reshape(x_c, (n_patch, n_c)), dims=1) ./ n_patch
    dr = abs.(p_r_old - p_r_new) ./ p_r_old
    dc = abs.(p_c_old - p_c_new) ./ p_c_old
    dr[isnan.(dr)] .= 0
    dc[isnan.(dc)] .= 0

	# return
	return dr, dc
end


function iterate_model_until_steady_state(x_state, x_r, x_c, n_r, n_c, M_inc, er, ec, cr, cc, n, n_patch, tmax, tmin, tol, k)

	# iterate until state
	for g=2:(tmax+1)

		# make copies of variables
		x_r_old = copy(x_r)
		x_c_old = copy(x_c)

		# update each patch sequentially
		x_r, x_c = update_patch_sequentially(M_inc, n, n_r, n_c, x_state, x_r, x_c, x_r_old, x_c_old, er, ec, cr, cc, k)

		# check convergence
		if g > tmin

			# check model convergence
			dr, dc = check_model_convergence(n_patch, n_r, n_c, x_r_old, x_c_old, x_r, x_c)
			if all(i -> i < tol, dr) && all(i -> i < tol, dc)
				break
			end
		end

	end # g loop

	# return
	return x_r, x_c
end
