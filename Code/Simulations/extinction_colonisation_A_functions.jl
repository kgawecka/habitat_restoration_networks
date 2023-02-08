# FUNCTIONS USED TO CALCULATE EXTINCTION AND COLONISATION PROBABILITIES - ANTAGONISM

function get_consumer_partners_in_current_patch(M_inc, x_c_old, i, j, h)

	# get consumer partners of resource h
	c_part = findall(M_inc[h,:].===1)

	# get consumers present in current patch
	c_current_all = findall(x_c_old[i,j,:].===1)

	# get consumer partners present in current patch
	c_part_sp = c_current_all[findall(x-> x in c_part, c_current_all)]

	# return
	return c_part_sp
end


function compute_resource_extinction_prob(er, c_part_sp)

	# initialise "not extinct" probability
	pe_r = 1-er

	# if at least one consumer partner present in current patch
	if length(c_part_sp) > 0

		i_count = 0

		for m=c_part_sp

			i_count = i_count+1
			prod = 1 - er / i_count
			pe_r = pe_r * prod

		end # m loop
	end

	# return
	return pe_r
end


function resource_extinction_prob(er, x_r, c_part_sp, i, j, h)

	# "not extinct" probability
	pe_r = compute_resource_extinction_prob(er, c_part_sp)

	if rand() > pe_r
		x_r[i,j,h] = 0
	end

	# return
	return x_r
end


function resource_extinctions(M_inc, er, x_r, x_c_old, i, j, h)

	# get resource partners in current patch
	c_part_sp = get_consumer_partners_in_current_patch(M_inc, x_c_old, i, j, h)

	# calculate extinction probability and update grid
	x_r = resource_extinction_prob(er, x_r, c_part_sp, i, j, h)

	# return
	return x_r
end


function get_resources_in_neighbouring_patches(x_r_old, neigh, h)

	# number of neighbouring patches with resource h
	n_neigh_r = 0

	# check neighbours
	for f = 1:size(neigh,1)

		i_n = neigh[f,1]
		j_n = neigh[f,2]

		# if resource h present in neighbouring patch
		if x_r_old[i_n, j_n, h] === 1
			n_neigh_r = n_neigh_r + 1
		end

	end # f loop

	# return
	return n_neigh_r
end


function resource_colonisations(cr, x_r_old, x_r, neigh, i, j, h)

	# number of neighbouring patches with resource h
	n_neigh_r = get_resources_in_neighbouring_patches(x_r_old, neigh, h)

	# "not colonised" probability
	pc_r = (1 - cr)^n_neigh_r

	if rand() > pc_r
		x_r[i,j,h] = 1
	end

	# return
	return x_r
end


function resource_extinctions_and_colonisations(M_inc, er, cr, n_r, x_r_old, x_r, x_c_old, neigh, i, j)

	# resource extinctions and colonisations
	for h=1:n_r

    	# if species present -> extinctions
    	if x_r_old[i, j, h] === 1
			x_r = resource_extinctions(M_inc, er, x_r, x_c_old, i, j, h)
        	continue
    	end

    	# if species absent -> colonisations
    	if x_r_old[i,j,h] === 0
			x_r = resource_colonisations(cr, x_r_old, x_r, neigh, i, j, h)
    	end
    end

	# return
	return x_r
end



function get_resource_partners_in_current_patch(M_inc, x_r_old, i, j, h)

	# get resource partners of consumer h
	r_part = findall(M_inc[:,h].===1)

	# get resources present in current patch
	r_current_all = findall(x_r_old[i,j,:].===1)

	# get resource partners present in current patch
	r_part_sp = r_current_all[findall(x-> x in r_part, r_current_all)]

	# return
	return r_part_sp
end


function compute_consumer_extinction_prob(ec, r_part_sp)

	# initialise extinction probability
	pe_c = ec

	# if at least one resource partner present in current patch
	if length(r_part_sp) > 0

		i_count = 0

		for m=r_part_sp

			i_count = i_count+1
			prod = 1 - ec / i_count
			pe_c = pe_c * prod

		end # m loop
	end

	# return
	return pe_c
end


function consumer_extinction_prob(ec, x_c, r_part_sp, i, j, h)

	# extinction probability
	pe_c = compute_consumer_extinction_prob(ec, r_part_sp)

	if rand() < pe_c
		x_c[i,j,h] = 0
	end

	# return
	return x_c
end


function consumer_extinctions(M_inc, ec, x_c, x_r_old, i, j, h)

	# get resource partners in current patch
	r_part_sp = get_resource_partners_in_current_patch(M_inc, x_r_old, i, j, h)

	# calculate extinction probability and update grid
	x_c = consumer_extinction_prob(ec, x_c, r_part_sp, i, j, h)

	# return
	return x_c
end


function get_consumers_in_neighbouring_patches(x_c_old, neigh, h)

	# number of neighbouring patches with consumer h
	n_neigh_c = 0

	# check neighbours
	for f = 1:size(neigh,1)

		i_n = neigh[f,1]
		j_n = neigh[f,2]

		# if consumer h present in neighbouring patch
		if x_c_old[i_n, j_n, h] === 1
			n_neigh_c = n_neigh_c + 1
		end

	end # f loop

	# return
	return n_neigh_c
end


function consumer_colonisations(cc, x_c_old, x_c, neigh, i, j, h)

	# number of neighbouring patches with consumer h
	n_neigh_c = get_consumers_in_neighbouring_patches(x_c_old, neigh, h)

	# "not colonised" probability
	pc_c = (1 - cc)^n_neigh_c

	if rand() > pc_c
		x_c[i,j,h] = 1
	end

	# return
	return x_c
end


function consumer_extinctions_and_colonisations(M_inc, ec, cc, n_c, x_c_old, x_c, x_r_old, neigh, i, j)

	# consumer extinctions and colonisations
	for h=1:n_c

    	# if species present -> extinctions
    	if x_c_old[i, j, h] === 1
			x_c = consumer_extinctions(M_inc, ec, x_c, x_r_old, i, j, h)
        	continue
    	end

    	# if species absent -> colonisations
    	if x_c_old[i,j,h] === 0
			x_c = consumer_colonisations(cc, x_c_old, x_c, neigh, i, j, h)
    	end
    end

	# return
	return x_c
end
