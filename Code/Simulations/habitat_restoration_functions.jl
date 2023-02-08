# FUNCTIONS USED TO SIMULATE HABITAT RESTORATION

function habitat_restoration_random(n, n_patch, dn, x_state, k)

	# create matrix with patch numbers (col1) and states (col2)
	m_state = [1:n_patch  vec(x_state[:,:,k-1])]

	# matrix with destroyed patches only
	m_destroyed = m_state[m_state[:,2].===0,1]

	# select dn patches at random and change their state to 1 (pristine)
	patch_d = sample(m_destroyed, dn, replace=false)
	m_state[in.(m_state[:,1], Ref(patch_d)),2] .= 1

	# new grid
	x_state[:,:,k] = reshape(m_state[:,2], (n,n))

	# return
	return x_state
end


function habitat_restoration_nonrandom(n, n_patch, dn, patch_neigh, x_state, x_r, x_c, k)

	# create matrix with patch numbers (col1) and states (col2)
	m_state = [1:n_patch  vec(x_state[:,:,k-1])]

	# matrix with destroyed patches only
	m_destroyed = m_state[m_state[:,2].===0,1]

	# grid with number of species in each patch
	x_species = sum(x_r, dims=3) + sum(x_c, dims=3)

	# create matrix with patch numbers (col1) and number of species (col2)
	m_species = [1:n_patch  vec(x_species)]

	# matrix with occupied patches only
	m_species = m_species[m_species[:,2].!==0,1]

	# neighbours of occupied patches
	neigh_occ = unique(patch_neigh[in.(patch_neigh.patch_no, (m_species,)),"neighbour"])

	# neighbours of occupied patches that are destroyed
    neigh_destroyed = neigh_occ[in.(neigh_occ, Ref(m_destroyed))]

	if length(neigh_destroyed)>dn
		# select dn patches at random and change their state to 1 (pristine)
		patch_r = sample(neigh_destroyed, dn, replace=false)
		m_state[in.(m_state[:,1], Ref(patch_r)),2] .= 1

	else
		# change state of all patches to 1 (pristine)
		m_state[in.(m_state[:,1], Ref(neigh_destroyed)),2] .= 1

		if length(neigh_destroyed)<dn

			# number of remaining "to be restored" patches
			nr = dn - length(neigh_destroyed)

			# matrix with destroyed patches only
			m_destroyed = m_state[m_state[:,2].===0,1]

			# select nr patches at random and change their state to 1 (pristine)
			patch_r = sample(m_destroyed, nr, replace=false)
			m_state[in.(m_state[:,1], Ref(patch_r)),2] .= 1

		end

	end

	# new grid
	x_state[:,:,k] = reshape(m_state[:,2], (n,n))

	# return
	return x_state
end


function habitat_restoration_target_diversity(n, n_patch, dn, patch_neigh, x_state, x_r, x_c, k)

	# create matrix with patch numbers (col1) and states (col2)
	m_state = [1:n_patch  vec(x_state[:,:,k-1])]

	# matrix with destroyed patches only
	m_destroyed = m_state[m_state[:,2].===0,1]

	# grid with number of species in each patch
	x_species = sum(x_r, dims=3) + sum(x_c, dims=3)

	# create data frame with patch numbers (col1) and number of species (col2)
	m_species = DataFrame(neighbour=1:n_patch, n_species=vec(x_species))

	# neighbours of destroyed patches
	neigh_destroyed = patch_neigh[in.(patch_neigh.patch_no, (m_destroyed,)),:]

	# combine with number of species in neighbouring patches
	neigh_destroyed = leftjoin(neigh_destroyed, m_species, on="neighbour")

	# get total number of species in all neighbouring patches
	neigh_destroyed = combine(groupby(neigh_destroyed, :patch_no), :n_species => sum)

	# sort by reducing number of species in neighbouring patches
	neigh_destroyed = neigh_destroyed[sortperm(neigh_destroyed[:,"n_species_sum"], rev=true),:]
	
	# patches with at least on species in neighbouring patches
	neigh_destroyed_species = neigh_destroyed[neigh_destroyed.n_species_sum.>0,"patch_no"]
	
	if length(neigh_destroyed_species)>dn
	  # select dn top patches change their state to 1 (pristine)
	  patch_r = neigh_destroyed[1:dn,"patch_no"]
	  m_state[in.(m_state[:,1], Ref(patch_r)),2] .= 1
	  
	else
	  # change state of all patches to 1 (pristine)
		m_state[in.(m_state[:,1], Ref(neigh_destroyed_species)),2] .= 1

		if length(neigh_destroyed_species)<dn

			# number of remaining "to be restored" patches
			nr = dn - length(neigh_destroyed_species)

			# matrix with destroyed patches only
			m_destroyed = m_state[m_state[:,2].===0,1]

			# select nr patches at random and change their state to 1 (pristine)
			patch_r = sample(m_destroyed, nr, replace=false)
			m_state[in.(m_state[:,1], Ref(patch_r)),2] .= 1

		end

	end

	# new grid
	x_state[:,:,k] = reshape(m_state[:,2], (n,n))

	# return
	return x_state
end


function habitat_restoration_target_species(n, n_patch, dn, patch_neigh, x_state, x_r, x_c, target_guild, target_species, k)

	# create matrix with patch numbers (col1) and states (col2)
	m_state = [1:n_patch  vec(x_state[:,:,k-1])]

	# matrix with destroyed patches only
	m_destroyed = m_state[m_state[:,2].===0,1]

	# create matrix with patch numbers (col1) and presence of target species (col2)
	if target_guild === "resource"
		m_species = [1:n_patch  vec(x_r[:,:,target_species])]
	elseif target_guild === "consumer"
		m_species = [1:n_patch  vec(x_c[:,:,target_species])]
	end

	# matrix with occupied patches only
	m_species = m_species[m_species[:,2].!==0,1]

	# neighbours of occupied patches
	neigh_occ = unique(patch_neigh[in.(patch_neigh.patch_no, (m_species,)),"neighbour"])

	# neighbours of occupied patches that are destroyed
    neigh_destroyed = neigh_occ[in.(neigh_occ, Ref(m_destroyed))]

	if length(neigh_destroyed)>dn
		# select dn patches at random and change their state to 1 (pristine)
		patch_r = sample(neigh_destroyed, dn, replace=false)
		m_state[in.(m_state[:,1], Ref(patch_r)),2] .= 1

	else
		# change state of all patches to 1 (pristine)
		m_state[in.(m_state[:,1], Ref(neigh_destroyed)),2] .= 1

		if length(neigh_destroyed)<dn # restore patches adjacent to any occupied patches

			# number of remaining "to be restored" patches
			nr = dn - length(neigh_destroyed)

			# matrix with destroyed patches only
			m_destroyed = m_state[m_state[:,2].===0,1]

			# grid with number of species in each patch
			x_species = sum(x_r, dims=3) + sum(x_c, dims=3)

			# create matrix with patch numbers (col1) and number of species (col2)
			m_species = [1:n_patch  vec(x_species)]

			# matrix with occupied patches only
			m_species = m_species[m_species[:,2].!==0,1]

			# neighbours of occupied patches
			neigh_occ = unique(patch_neigh[in.(patch_neigh.patch_no, (m_species,)),"neighbour"])

			# neighbours of occupied patches that are destroyed
		    neigh_destroyed = neigh_occ[in.(neigh_occ, Ref(m_destroyed))]

			if length(neigh_destroyed)>nr
				# select nr patches at random and change their state to 1 (pristine)
				patch_r = sample(m_destroyed, nr, replace=false)
				m_state[in.(m_state[:,1], Ref(patch_r)),2] .= 1

			else
				# change state of all patches to 1 (pristine)
				m_state[in.(m_state[:,1], Ref(neigh_destroyed)),2] .= 1

				if length(neigh_destroyed)<nr # restore patches at random

					# number of remaining "to be restored" patches
					nr = nr - length(neigh_destroyed)

					# matrix with destroyed patches only
					m_destroyed = m_state[m_state[:,2].===0,1]

					# select nr patches at random and change their state to 1 (pristine)
					patch_r = sample(m_destroyed, nr, replace=false)
					m_state[in.(m_state[:,1], Ref(patch_r)),2] .= 1

				end

			end

		end

	end

	# new grid
	x_state[:,:,k] = reshape(m_state[:,2], (n,n))

	# return
	return x_state
end