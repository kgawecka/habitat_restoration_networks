# FUNCTIONS USED TO CALCULATE EXTINCTION AND COLONISATION PROBABILITIES - MUTUALISM

function resource_extinctions(er, x_r, i, j, h)

	if rand() < er
		x_r[i,j,h] = 0
	end

	# return grid of plants and grid of plant trait values
	return x_r
end


function get_consumer_partners_in_current_patch(M_inc, h)

	# find consumer partners of resource h
	c_part = findall(M_inc[h,:].===1)

	# return
	return c_part
end


function compute_resource_colonisation_prob(pc_r, cr, c_part_sp)

	# "not colonised" probability
	j_count = 0

	for l=c_part_sp

		j_count = j_count+1
		prod = 1 - cr / j_count
		pc_r = pc_r * prod

	end # l loop

	# return
	return pc_r
end


function resource_colonisation_prob(cr, x_r_old, x_r, x_c_old, c_part, neigh, i, j, h)

	# initialise "not colonised" probability
	pc_r = 1

	# check neighbours
    for f=1:size(neigh,1)

		# neighbour coordinates
    	i_n = neigh[f,1]
    	j_n = neigh[f,2]

		# if resource h absent from neighbouring patch
    	if x_r_old[i_n,j_n,h] === 0
    		continue

    	# if resource h present in neighbouring patch
		else

			# consumers present in neighbouring patch
			c_neigh_all = findall(x_c_old[i_n,j_n,:].===1)

			# consumer partners present in neighbouring patch
			c_part_sp = c_neigh_all[findall(x-> x in c_part, c_neigh_all)]

			# "not colonised" probability
			if length(c_part_sp) !== 0 # if resource h and its partner(s) present in neighbouring patch
				pc_r = compute_resource_colonisation_prob(pc_r, cr, c_part_sp)
			end

    	end

    end # f loop

	if rand() > pc_r
		x_r[i,j,h] = 1
	end

	# return
	return(x_r)
end


function resource_colonisations(M_inc, cr, x_r_old, x_r, x_c_old, neigh, i, j, h)

	# find consumer partners of resource
	c_part = get_consumer_partners_in_current_patch(M_inc, h)

	x_r = resource_colonisation_prob(cr, x_r_old, x_r, x_c_old, c_part, neigh, i, j, h)

	# return
	return x_r
end


function resource_extinctions_and_colonisations(M_inc, er, cr, n_r, x_r_old, x_r, x_c_old, neigh, i, j)

	# resource extinctions and colonisations
    for h=1:n_r

    	# if species present -> extinctions
        if x_r_old[i,j,h] === 1
    		x_r = resource_extinctions(er, x_r, i, j, h)
        continue
        end

        # if species absent -> colonisations
        if x_r_old[i,j,h] === 0
        	x_r = resource_colonisations(M_inc, cr, x_r_old, x_r, x_c_old, neigh, i, j, h)
        end

    end

	# return
	return x_r
end



function consumer_extinctions(ec, x_c, i, j, h)

	if rand() < ec
		x_c[i,j,h] = 0
	end

	# return
	return x_c
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


function compute_consumer_colonisation_prob(pc_c, cc, r_part_sp)

	# "not colonised" probability
	i_count = 0

	for m=r_part_sp

		i_count = i_count+1
		prod = 1 - cc / i_count
		pc_c = pc_c * prod

	end # m loop

	# return
	return pc_c
end


function consumer_colonisation_prob(cc, x_c_old, x_c, r_part_sp, neigh, i, j, h)

	# initialise "not colonised" probability
	pc_c = 1

	for f = 1:size(neigh,1)

		i_n = neigh[f,1]
		j_n = neigh[f,2]

		# if consumer h absent from neighbouring patch
		if x_c_old[i_n, j_n, h] === 0
			continue

		# if consumer h present in neighbouring patch
		else
			# "not colonised" probability
			pc_c = compute_consumer_colonisation_prob(pc_c, cc, r_part_sp)
		end

	end # f loop

	if rand() > pc_c
		x_c[i,j,h] = 1
	end

	# return
	return x_c
end


function consumer_colonisations(M_inc, cc, x_c_old, x_c, x_r_old, neigh, i, j, h)

	# get resource partners in current patch
	r_part_sp = get_resource_partners_in_current_patch(M_inc, x_r_old, i, j, h)

	# check neighbours
	# if at least one resource partner present in current patch
	if length(r_part_sp) > 0
		x_c = consumer_colonisation_prob(cc, x_c_old, x_c, r_part_sp, neigh, i, j, h)
	end

	# return
	return x_c
end


function consumer_extinctions_and_colonisations(M_inc, ec, cc, n_c, x_c_old, x_c, x_r_old, neigh, i, j)

	# consumer extinctions and colonisations
	for h=1:n_c

    	# if species present -> extinctions
    	if x_c_old[i, j, h] === 1
			x_c = consumer_extinctions(ec, x_c, i, j, h)
        	continue
    	end

    	# if species absent -> colonisations
    	if x_c_old[i,j,h] === 0
			x_c = consumer_colonisations(M_inc, cc, x_c_old, x_c, x_r_old, neigh, i, j, h)
    	end
    end

	# return
	return x_c
end
