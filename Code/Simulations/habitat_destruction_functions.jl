# FUNCTIONS USED TO SIMULATE HABITAT DESTRUCTION

function habitat_destruction(n, n_patch, dn, x_state, x_r, x_c, k)

    # create matrix with patch numbers (col1) and states (col2)
    m_state = [1:n_patch  vec(x_state[:,:,k])]

    # matrix with pristine patches only
    m_pristine = m_state[m_state[:,2].!==0,:]

    # select dn patches at random and change their state to 0 (destroyed)
    patch_d = sample(m_pristine[:,1], dn, replace=false)
    m_state[in.(m_state[:,1], Ref(patch_d)),2] .= 0

    # new grid
    x_state[:,:,k+1] = reshape(m_state[:,2], (n,n))

    # species in destroyed patches become extinct
    x_r = x_r .* x_state[:,:,k+1]
    x_c = x_c .* x_state[:,:,k+1]

	# return
	return x_state, x_r, x_c

end


function habitat_destruction_nonrandom(n, n_patch, dn, x_state, x_r, x_c, k, patch_neigh)

    # create matrix with patch numbers (col1) and states (col2)
    m_state = [1:n_patch  vec(x_state[:,:,k])]

    # matrix with pristine patches only
    m_pristine = m_state[m_state[:,2].!==0,:]

    if k===1  # first destruction step - destroy 100 patches at random first
        # select 100 patches at random and change their state to 0 (destroyed)
        patch_d = sample(m_pristine[:,1], 100, replace=false)
        m_state[in.(m_state[:,1], Ref(patch_d)),2] .= 0

        # number of remaining "to be destroyed" patches
        n_d = dn - 100

    else     # other destruction steps - destroy dn patches
        n_d = dn
    end

    # matrix with destroyed patches only
    m_destroyed = m_state[m_state[:,2].===0,:]

    # neighbours of destroyed patches
    destroyed_neigh = unique(patch_neigh[in.(patch_neigh.patch_no, (m_destroyed[:,1],)),"neighbour"])

    # matrix with pristine patches only
    m_pristine = m_state[m_state[:,2].!==0,:]

    # neighbours of destroyed patches that are pristine
    neigh_pristine = m_pristine[in.(m_pristine[:,1], (destroyed_neigh,)),1]

    if length(neigh_pristine)>n_d
        # select n_d patches at random and change their state to 0 (destroyed)
        patch_d = sample(neigh_pristine, n_d, replace=false)
        m_state[in.(m_state[:,1], Ref(patch_d)),2] .= 0

    else
        # select all patches and change their state to 0 (destroyed)
        m_state[in.(m_state[:,1], Ref(neigh_pristine)),2] .= 0

        if length(neigh_pristine)<n_d
            # number of remaining "to be destroyed" patches
            n_d = n_d - length(neigh_pristine)

            # matrix with destroyed patches only
            m_destroyed = m_state[m_state[:,2].===0,:]

            # neighbours of destroyed patches
            destroyed_neigh = unique(patch_neigh[in.(patch_neigh.patch_no, (m_destroyed[:,1],)),"neighbour"])

            # matrix with pristine patches only
            m_pristine = m_state[m_state[:,2].!==0,:]

            # neighbours of destroyed patches that are pristine
            neigh_pristine = m_pristine[in.(m_pristine[:,1], (destroyed_neigh,)),1]

            if length(neigh_pristine)>n_d
                # select n_d patches at random and change their state to 0 (destroyed)
                patch_d = sample(neigh_pristine, n_d, replace=false)
                m_state[in.(m_state[:,1], Ref(patch_d)),2] .= 0

            else
                # select all patches and change their state to 0 (destroyed)
                m_state[in.(m_state[:,1], Ref(neigh_pristine)),2] .= 0

                if length(neigh_pristine)<n_d
                    # number of remaining "to be destroyed" patches
                    n_d = n_d - length(neigh_pristine)

                    # matrix with pristine patches only
                    m_pristine = m_state[m_state[:,2].!==0,:]

                    # select n_d patches at random and change their state to 0 (destroyed)
                    patch_d = sample(m_pristine[:,1], n_d, replace=false)
                    m_state[in.(m_state[:,1], Ref(patch_d)),2] .= 0
                end
            end
        end
    end
    
    # new grid
    x_state[:,:,k+1] = reshape(m_state[:,2], (n,n))

    # species in destroyed patches become extinct
    x_r = x_r .* x_state[:,:,k+1]
    x_c = x_c .* x_state[:,:,k+1]

	# return
	return x_state, x_r, x_c
	
end