include("src/QSN.jl")

using .QSN
using SparseArrays
using LinearAlgebra
using ITensors
using Plots
using GraphPlot
using GraphRecipes

mutable struct Graph
    vertices::Array{Int64}
    edges::Dict{Tuple{Int64, Int64}, Union{Float64, Array{Float64, 1}}}

    Graph() = new([], Dict())

    Graph(n::Int64) = new([i for i in 1:n], Dict())
end

add_vertex!(g::Graph) = push!(g.vertices, length(g.vertices)+1)
add_vertices!(g::Graph, n::Int64) = append!(g.vertices, collect(length(g.vertices)+1:length(g.vertices)+1+n))

function add_edge!(g::Graph, v1::Int64, v2::Int64, w::Float64)
    if haskey(g.edges, (v1, v2))
        if typeof(g.edges[(v1, v2)]) <: Vector
            push!(g.edges[(v1, v2)], w)
        else
            g.edges[(v1, v2)] = Float64[g.edges[(v1, v2)], w]
        end
    else
        g.edges[(v1, v2)] = w
    end
end

function rem_edge!(g::Graph, v1::Int64, v2::Int64, w::Float64)
    if typeof(g.edges[(v1, v2)]) <: Vector
        deleteat!(g.edges[(v1, v2)], findfirst(x -> x == w, g.edges[(v1, v2)]))
        if length(g.edges[(v1, v2)]) == 1
            g.edges[(v1, v2)] = g.edges[(v1, v2)][1]
        end
    else
        delete!(g.edges, (v1, v2))
    end
end

function get_weights(g::Graph, v1::Int64, v2::Int64)
    return g.edges[(v1, v2)]
end

function edges(g::Graph)
    es = []
    for edge in keys(g.edges)
        if typeof(g.edges[edge]) <: Vector
            append!(es, [edge for w in g.edges[edge]])
        else
            push!(es, edge)
        end
    end
    return es
end

function adjacency_list_and_labels(g::Graph)
    adj_list = [Int[] for _ in 1:length(g.vertices)]
    edge_labels = Dict{Tuple{Int, Int}, String}()
    
    for (v1, v2) in keys(g.edges)
        weights = get_weights(g, v1, v2)
        if typeof(weights) <: Vector
            for w in weights
                append!(adj_list[v1], v2)
                edge_labels[(v1, v2)] = string(weights)
            end
        else
            append!(adj_list[v1], v2)
            edge_labels[(v1, v2)] = string(weights)
        end
    end

    return adj_list, edge_labels
end

# Plot the graph using graphplot
function plot_graph(g::Graph)
    adj_list, edge_labels = adjacency_list_and_labels(g)
    names = string.(1:length(g.vertices))
    # Calculate grid layout dimensions
    ncols = ceil(Int, sqrt(length(g.vertices)))
    nrows = ceil(Int, length(g.vertices) / ncols)
    
    # Generate coordinates for grid layout
    coords_x = repeat(1:ncols, inner=nrows)
    coords_y = repeat(1:nrows, outer=ncols)[1:length(g.vertices)]

    coords_x = map(Float64, coords_x)
    coords_y = map(Float64, coords_y)
    
    # Plot the graph with specified parameters
    graphplot(adj_list, 
              names=names, 
              nodeshape=:circle, 
              self_edge_size=0.4,
              nodesize = 0.2,
              arrow=false,
              edgelabel=Dict(((u, v), edge_labels[(u, v)]) for (u, v) in keys(edge_labels)),
              coords=(coords_x, coords_y)
    )
end


# """
#     initial_state(::Array{Int64,1})

# Returns an initial state in MPS form for a given qubit order

# Parameters:

#     qubits::Array{Int64,1}: an array of 0s and 1s, describing the state

# Output:

#     state::Array{ITensor,1}: array of tensors in the state
# """
# function initial_state(qubits::Array{Int64,1})

#     # Initializing number of qubits
#     N = length(qubits)

#     # Generate MPS with N sites and site dimension 2
#     state = MPS(siteinds(2, N))

#     # Site 1
#     if qubits[1] == 0
#         state[1] = N == 1 ? ITensor([1.0+0.0im;0.0+0.0im], settags(inds(state[1])[1], "qubit=1,layer=0")) : 
#                             ITensor([1.0+0.0im;0.0+0.0im;;], (settags(inds(state[1])[2], "qubit=1,layer=0"), 
#                                                                 settags(inds(state[1])[1], "Link,n=1")))
#     elseif qubits[1] == 1
#         state[1] = N == 1 ? ITensor([0.0+0.0im;1.0+0.0im], settags(inds(state[1])[1], "qubit=1,layer=0")) : 
#                             ITensor([0.0+0.0im;1.0+0.0im;;], (settags(inds(state[1])[2], "qubit=1,layer=0"), 
#                                                                 settags(inds(state[1])[1], "Link,n=1")))
#     end

#     # Sites from 2:(N-1)
#     for i in 2:(N-1)
#         if qubits[i] == 0
#             state[i] = ITensor([1.0+0.0im 0.0+0.0im;;;], (settags(inds(state[i])[1], "Link,n=$(i-1)"), settags(inds(state[i])[2], "layer=0,qubit=$i"), settags(inds(state[i])[3], "Link,n=$(i)")))
#         elseif qubits[i] == 1
#             state[i] = ITensor([0.0+0.0im 1.0+0.0im;;;], (settags(inds(state[i])[1], "Link,n=$(i-1)"), settags(inds(state[i])[2], "layer=0,qubit=$i"), settags(inds(state[i])[3], "Link,n=$(i)")))
#         end
#     end

#     # Site N
#     if qubits[end] == 0 && N != 1
#         state[end] = ITensor([1.0+0.0im 0.0+0.0im], (settags(inds(state[end])[1], "Link,n=$(N-1)"), settags(inds(state[end])[2], "layer=0,qubit=$N")))
#     elseif qubits[end] == 1 && N != 1
#         state[end] = ITensor([0.0+0.0im 1.0+0.0im], (settags(inds(state[end])[1], "Link,n=$(N-1)"), settags(inds(state[end])[2], "layer=0,qubit=$N")))
#     end

#     return state
# end

# """
#     ket0(::Int64)

# Returns an N-qubit |0...0> state
# """
# ket0(N::Int64) = initial_state([0 for i in 1:N])


function entanglement_swapping!(g, v1, v2, v3)
    α = get_weights(g, v1, v2)
    β = get_weights(g, v2, v3)

    λ = 0.5*(1 + √(1 - 16*α*(1-α)*β*(1-β)))

    rem_edge!(g, v1, v2, α)
    rem_edge!(g, v2, v3, β)
    add_edge!(g, v1, v3, λ)
end

function entanglement_distillation!(g, v1, v2)
    α, β = get_weights(g, v1, v2)[1:2]

    λ = max(0.5, α*β)

    rem_edge!(g, v1, v2, α)
    rem_edge!(g, v1, v2, β)
    add_edge!(g, v1, v2, λ)
end

g = Graph()

add_vertices!(g, 16)

add_edge!(g, 1, 2, rand(0.5:0.1:0.9))
add_edge!(g, 2, 3, rand(0.5:0.1:0.9))
add_edge!(g, 3, 4, rand(0.5:0.1:0.9))
add_edge!(g, 1, 5, rand(0.5:0.1:0.9))
add_edge!(g, 2, 6, rand(0.5:0.1:0.9))
add_edge!(g, 3, 7, rand(0.5:0.1:0.9))
add_edge!(g, 4, 8, rand(0.5:0.1:0.9))
add_edge!(g, 5, 6, rand(0.5:0.1:0.9))
add_edge!(g, 6, 7, rand(0.5:0.1:0.9))
add_edge!(g, 7, 8, rand(0.5:0.1:0.9))
add_edge!(g, 5, 9, rand(0.5:0.1:0.9))
add_edge!(g, 6, 10, rand(0.5:0.1:0.9))
add_edge!(g, 7, 11, rand(0.5:0.1:0.9))
add_edge!(g, 8, 12, rand(0.5:0.1:0.9))
add_edge!(g, 9, 10, rand(0.5:0.1:0.9))
add_edge!(g, 10, 11, rand(0.5:0.1:0.9))
add_edge!(g, 11, 12, rand(0.5:0.1:0.9))
add_edge!(g, 9, 13, rand(0.5:0.1:0.9))
add_edge!(g, 10, 14, rand(0.5:0.1:0.9))
add_edge!(g, 11, 15, rand(0.5:0.1:0.9))
add_edge!(g, 12, 16, rand(0.5:0.1:0.9))
add_edge!(g, 13, 14, rand(0.5:0.1:0.9))
add_edge!(g, 14, 15, rand(0.5:0.1:0.9))
add_edge!(g, 15, 16, rand(0.5:0.1:0.9))

entanglement_swapping!(g, 1, 2, 6)
entanglement_swapping!(g, 1, 5, 6)
entanglement_distillation!(g, 1, 6)

entanglement_swapping!(g, 6, 7, 11)
entanglement_swapping!(g, 6, 10, 11)
entanglement_distillation!(g, 6, 11)

entanglement_swapping!(g, 11, 12, 16)
entanglement_swapping!(g, 11, 15, 16)
entanglement_distillation!(g, 11, 16)

plt = plot_graph(g)
savefig(plt, "graph.png")