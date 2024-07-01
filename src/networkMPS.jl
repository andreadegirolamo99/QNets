function _balanced_sort(arr::Array{Float64, 1})
    nodes = sortperm(arr)
    arr = arr[nodes]

    res = zeros(length(arr))
    nodes_perm = zeros(Int64, length(arr))
    count_up = 1
    count_down = length(arr)
    for i in eachindex(res)
        if i % 2 == 1
            res[count_up] = arr[i]
            nodes_perm[count_up] = nodes[i]
            count_up += 1
        else
            res[count_down] = arr[i]
            nodes_perm[count_down] = nodes[i]
            count_down -= 1
        end
    end
    
    return nodes_perm
end

function initial_state(qubits::Array{Int64,1}, photons::Array{Photon,1}, p_to_n::Dict{UInt64, UInt64})
    # Initializing number of qubits
    N = length(qubits)

    # Generate MPS with N sites and site dimension 2
    state = MPS(siteinds(2, N))

    # Site 1
    if qubits[1] == 0
        state[1] = N == 1 ? ITensor([1.0+0.0im;0.0+0.0im], 
                                settags(inds(state[1])[1], "Node=$(p_to_n[photons[1].id]),Photon=$(photons[1].id % 1000)")) : 
                            ITensor([1.0+0.0im;0.0+0.0im;;], 
                                (settags(inds(state[1])[2], "Node=$(p_to_n[photons[1].id]),Photon=$(photons[1].id % 1000)"), 
                                settags(inds(state[1])[1], "Link,n=1")))
    elseif qubits[1] == 1
        state[1] = N == 1 ? ITensor([0.0+0.0im;1.0+0.0im], 
                                settags(inds(state[1])[1], "Node=$(p_to_n[photons[1].id]),Photon=$(photons[1].id % 1000)")) : 
                            ITensor([0.0+0.0im;1.0+0.0im;;], 
                                (settags(inds(state[1])[2], "Node=$(p_to_n[photons[1].id]),Photon=$(photons[1].id % 1000)"), 
                                settags(inds(state[1])[1], "Link,n=1")))
    end

    # Sites from 2:(N-1)
    for i in 2:(N-1)
        if qubits[i] == 0
            state[i] = ITensor([1.0+0.0im 0.0+0.0im;;;], 
                        (settags(inds(state[i])[1], "Link,n=$(i-1)"), 
                        settags(inds(state[i])[2], "Node=$(p_to_n[photons[i].id]),Photon=$(photons[i].id % 1000)"), 
                        settags(inds(state[i])[3], "Link,n=$(i)")))
        elseif qubits[i] == 1
            state[i] = ITensor([0.0+0.0im 1.0+0.0im;;;], 
                        (settags(inds(state[i])[1], "Link,n=$(i-1)"), 
                        settags(inds(state[i])[2], "Node=$(p_to_n[photons[i].id]),Photon=$(photons[i].id % 1000)"), 
                        settags(inds(state[i])[3], "Link,n=$(i)")))
        end
    end

    # Site N
    if qubits[end] == 0 && N != 1
        state[end] = ITensor([1.0+0.0im 0.0+0.0im], 
                        (settags(inds(state[end])[1], "Link,n=$(N-1)"), 
                        settags(inds(state[end])[2], "Node=$(p_to_n[photons[end].id]),Photon=$(photons[end].id % 1000)")))
    elseif qubits[end] == 1 && N != 1
        state[end] = ITensor([0.0+0.0im 1.0+0.0im], 
                        (settags(inds(state[end])[1], "Link,n=$(N-1)"), 
                        settags(inds(state[end])[2], "Node=$(p_to_n[photons[end].id]),Photon=$(photons[end].id % 1000)")))
    end

    return state
end

ket0(N::Int64, photons::Array{Photon,1}, p_to_n::Dict{UInt64, UInt64}) = initial_state([0 for i in 1:N], photons, p_to_n)

function connectMPS!(state::MPS, p1::Photon, p2::Photon, p_to_MPS::Dict)
    apply_unitary!(state, H, p_to_MPS[p1.id])
    println(state)
    apply_unitary!(state, CNOT, p_to_MPS[p1.id], p_to_MPS[p2.id])
end

function generate_network_MPS(qn::QuantumNetwork)
    node_order = _balanced_sort(closeness_centrality(Graph(adjMat(qn))))
    ordered_photons = [qn.nodes[i].photons[j] for i in node_order for j in eachindex(qn.nodes[i].photons)]

    p_to_n = Dict{UInt64,UInt64}()
    for n in qn.nodes
        for p in n.photons
            p_to_n[p.id] = n.id
        end
    end

    state = ket0(sum(n_photons(node) for node in qn.nodes), ordered_photons, p_to_n)
    p_to_MPS = Dict()
    for i in eachindex(ordered_photons)
        p_to_MPS[ordered_photons[i].id] = i
    end

    photon_connections = [connected_photons(qn.nodes[n1], qn.nodes[n2]) for (n1, n2) in qn.topology]
    for (p1,p2) in photon_connections
        connectMPS!(state, p1, p2, p_to_MPS)
    end

    return state
end