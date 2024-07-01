const I2 = [1.0+0.0im 0.0+0.0im; 
            0.0+0.0im 1.0+0.0im]                    # Identity gate
const X = [0.0+0.0im 1.0+0.0im; 
            1.0+0.0im 0.0+0.0im]                     # Pauli-X gate
const H = sqrt(0.5).*[1.0+0.0im 1.0+0.0im; 
                1.0+0.0im -1.0+0.0im]         # Hadamard gate
const SWAP = reshape([1.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im; 
                0.0+0.0im 0.0+0.0im 1.0+0.0im 0.0+0.0im; 
                0.0+0.0im 1.0+0.0im 0.0+0.0im 0.0+0.0im; 
                0.0+0.0im 0.0+0.0im 0.0+0.0im 1.0+0.0im], (2,2,2,2))       # SWAP gate

const CNOT = reshape([I2 ; zeros(2,2) ;; 
                zeros(2,2) ; X], (2,2,2,2))                                # CNOT gate

function apply_unitary!(state::MPS, U, sites...)
    if length(size(U)) == 2
        site = sites[1]
        prev_idx = inds(state[site])[site == 1 ? 1 : 2]
        new_idx = Index(2, tags(prev_idx))

        orthogonalize!(state, site)
        for i in 2:length(state)-1
            state[i] = ITensor(convert(Array, state[i].tensor), inds(state[i])[3], inds(state[i])[1], inds(state[i])[2])
        end
        state[end] = ITensor(convert(Array, state[end].tensor), inds(state[end])[2], inds(state[end])[1])

        state[site] *= ITensor(U, prev_idx, new_idx)
    elseif length(size(U)) == 4
        prev_idx_1 = inds(state[sites[1]])[site == 1 ? 1 : 2]
        prev_idx_2 = inds(state[sites[2]])[site == 1 ? 1 : 2]
        new_idx_1 = Index(2, tags(prev_idx_1))
        new_idx_2 = Index(2, tags(prev_idx_2))

        site1 = sites[1]
        if sites[2] > sites[1] + 1
            for site in sites[1]:sites[2]-2
                prev_idx_1 = inds(state[site])[site == 1 ? 1 : 2]
                prev_idx_2 = inds(state[site+1])[site == 1 ? 1 : 2]
                new_idx_1 = Index(2, tags(prev_idx_1))
                new_idx_2 = Index(2, tags(prev_idx_2))

                orthogonalize!(state, site)
                for i in 2:length(state)-1
                    state[i] = ITensor(convert(Array, state[i].tensor), inds(state[i])[3], inds(state[i])[1], inds(state[i])[2])
                end
                state[end] = ITensor(convert(Array, state[end].tensor), inds(state[end])[2], inds(state[end])[1])

                state[site] *= ITensor(SWAP, prev_idx_1, prev_idx_2, new_idx_1, new_idx_2)
            end
        end

        orthogonalize!(state, site1)
        

        if sites[2] > sites[1] + 1

        end
    end
    return nothing
end