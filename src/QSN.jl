module QSN

###################################################
# Exports
#
export
    # Photon
    Photon,
    isFree,
    # Node
    Node,
    n_photons,
    connections,
    free_photons,
    add_photon!,
    add_photons!,
    # QuantumNetwork
    QuantumNetwork,
    add_node!,
    get_node,
    connect!,
    isConnected,
    adjMat,
    # Network-MPS
    generate_network_MPS
#
###################################################


###################################################
# External packages
#
using ITensors
using Random
using SparseArrays
using Graphs
using StatsBase
using LinearAlgebra
using Distributed
#
###################################################


###################################################
# Internal files
#
include("quantumnetwork.jl")
include("networkMPS.jl")
include("operators.jl")
#
###################################################

end # module QSN