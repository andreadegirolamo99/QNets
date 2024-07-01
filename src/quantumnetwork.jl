const PHOTON_IDS = UInt64[]
function gen_photon_id()
    pid = rand(MersenneTwister(), UInt64)
    while pid in PHOTON_IDS
        pid = rand(MersenneTwister(), UInt64)
    end
    push!(PHOTON_IDS, pid)
    return pid
end

mutable struct Photon
    id::UInt64
    nextP::Union{Nothing, Photon}

    Photon() = new(gen_photon_id(), nothing)
end

isFree(p::Photon) = isnothing(p.nextP)

function Base.show(io::IO, p::Photon)
    print(io, "Photon(id=$(p.id % 1000)|$(isFree(p) ? "free" : "connected -> $(p.nextP.id % 1000)"))")
end

const NODE_IDS = UInt64[]
function gen_node_id()
    nid = isempty(NODE_IDS) ? 1 : NODE_IDS[end]+1
    push!(NODE_IDS, nid)
    return nid
end

mutable struct Node
    id::UInt64
    photons::Array{Photon, 1}

    Node() = new(gen_node_id(), [])
    Node(photons::Array{Photon, 1}) = new(gen_node_id(), photons)
    Node(N::Int64) = new(gen_node_id(), [Photon() for i in 1:N])
end

add_photon!(n::Node) = push!(n.photons, Photon())
add_photons!(n::Node, p::Int64) = append!(n.photons, [Photon() for i in 1:p])

n_photons(n::Node) = length(n.photons)
connections(n::Node) = count(p -> !isFree(p), n.photons)
free_photons(n::Node) = count(p -> isFree(p), n.photons)

hasphoton(n::Node, p::Photon) = p in n.photons

function Base.show(io::IO, n::Node)
    pidstr = join(["$(p.id % 1000)," for p in n.photons])[1:end-1]
    print(io, "Node(id=$(n.id)|photons=($pidstr))")
end

mutable struct QuantumNetwork
    nodes::Array{Node, 1}
    topology::Array{Tuple{UInt64, UInt64}, 1}

    QuantumNetwork() = new([], [])
    QuantumNetwork(nodes::Array{Node, 1}) = new(nodes, [])
end

node(photon::Photon, qn::QuantumNetwork) = first(filter(n -> hasphoton(n, photon), qn.nodes))

add_node!(qn::QuantumNetwork, node::Node) = push!(qn.nodes, node)
add_nodes!(qn::QuantumNetwork, nodes::Node) = append!(qn.nodes, nodes)

get_node(qn::QuantumNetwork, nid::UInt64) = first(filter(n -> n.id == nid, qn.nodes))

_no_free_photon_error(nid::UInt64) = error("node n.$nid has no free photons. Add a new photon before attempting a connection.")
_node_not_in_qn_error(nid::UInt64) = error("node n.$nid is not in the given network.")

function connect!(qn::QuantumNetwork, node1::Node, node2::Node)
    node1 in qn.nodes || _node_not_in_qn_error(node1.id)
    node2 in qn.nodes || _node_not_in_qn_error(node2.id)

    reduce(|, [isFree(p) for p in node1.photons]) || _no_free_photon_error(node1.id)
    reduce(|, [isFree(p) for p in node2.photons]) || _no_free_photon_error(node2.id)

    push!(qn.topology, (node1.id, node2.id))
    p1 = findfirst(p -> isFree(p), node1.photons)
    p2 = findfirst(p -> isFree(p) && p ≠ node1.photons[p1], node2.photons)
    node1.photons[p1].nextP = node2.photons[p2]
    node2.photons[p2].nextP = node1.photons[p1]

    return nothing
end

isConnected(n1::Node, n2::Node) = reduce(|, [p.nextP in n2.photons for p in n1.photons])

photons(qn::QuantumNetwork) = [n.photons[i] for n in qn.nodes for i in eachindex(n.photons)]

_n1n2_not_connected_error(n1::UInt64, n2::UInt64) = error("nodes $n1 and $n2 are not connected.")

function connected_photons(n1::Node, n2::Node)
    isConnected(n1, n2) || _n1n2_not_connected_error(n1.id, n2.id)
    return (first(filter(p -> p.nextP in n2.photons, n1.photons)), first(filter(p -> p.nextP in n1.photons, n2.photons)))
end

function adjMat(qn::QuantumNetwork)
    adj = spzeros(Int64, length(qn.nodes), length(qn.nodes))
    for i in eachindex(qn.nodes), j in eachindex(qn.nodes)
        if isConnected(qn.nodes[i], qn.nodes[j])
            adj[i,j] = 1
        end
    end
    return adj
end

function QuantumNetwork(adj::SparseMatrixCSC)
    nodes = [Node() for i in 1:size(adj)[1]]
    qn = QuantumNetwork(nodes)

    for i in eachindex(qn.nodes), j in i:length(qn.nodes)
        if adj[i, j] ≠ 0
            add_photon!(qn.nodes[i])
            add_photon!(qn.nodes[j])
            connect!(qn, qn.nodes[i], qn.nodes[j])
        end
    end

    return qn
end