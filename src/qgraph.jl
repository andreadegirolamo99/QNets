using DataStructures
using GraphPlot
using GraphRecipes
using Plots
using StatsPlots

mutable struct QGraph
    vertices::Array{Int64}
    edges::Dict{Tuple{Int64, Int64}, Union{Float64, Array{Float64, 1}}}

    QGraph() = new([], Dict())

    QGraph(n::Int64) = new([i for i in 1:n], Dict())
end

add_vertex!(g::QGraph) = push!(g.vertices, length(g.vertices)+1)
add_vertices!(g::QGraph, n::Int64) = append!(g.vertices, collect(length(g.vertices)+1:length(g.vertices)+1+n))

vertices(g::QGraph) = g.vertices
n_vertices(g::QGraph) = length(g.vertices)

function add_edge!(g::QGraph, v1::Int64, v2::Int64, w::Float64)
    v1, v2 = min(v1, v2), max(v1, v2)
    
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

function rem_edge!(g::QGraph, v1::Int64, v2::Int64, w::Float64)
    v1, v2 = min(v1, v2), max(v1, v2)
    if typeof(g.edges[(v1, v2)]) <: Vector
        deleteat!(g.edges[(v1, v2)], findfirst(x -> x == w, g.edges[(v1, v2)]))
        if length(g.edges[(v1, v2)]) == 1
            g.edges[(v1, v2)] = g.edges[(v1, v2)][1]
        end
    else
        delete!(g.edges, (v1, v2))
    end
end

get_weights(g::QGraph, v1::Int64, v2::Int64) = are_neighbors(g, v1, v2) ? (v1 < v2 ? g.edges[(v1, v2)] : g.edges[(v2, v1)]) : []
    
function edges(g::QGraph)
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

edges(g::QGraph, v::Int64) = filter(edge -> v in edge, keys(g.edges))

function get_neighbors(g::QGraph, v::Int64)
    neighbors = Int[]
    for (v1, v2) in keys(g.edges)
        if v1 == v
            push!(neighbors, v2)
        elseif v2 == v
            push!(neighbors, v1)
        end
    end
    return neighbors
end

are_neighbors(g::QGraph, v1::Int64, v2::Int64) = v2 in get_neighbors(g, v1)

function adjacency_list_and_labels(g::QGraph)
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

function adjMat(g::QGraph)
    n = length(g.vertices)
    rows = Int[]
    cols = Int[]
    vals = Float64[]
    
    for (v1, v2) in keys(g.edges)
        weights = get_weights(g, v1, v2)
        if typeof(weights) <: Vector
            for w in weights
                push!(rows, v1)
                push!(cols, v2)
                push!(vals, w)
            end
        else
            push!(rows, v1)
            push!(cols, v2)
            push!(vals, weights)
        end
    end

    return sparse(rows, cols, vals, n, n)
end

# Plot the graph using graphplot
function plot_graph(g::QGraph)
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
              coords=(coords_x, coords_y)
    )
end

function distance(g::QGraph, v::Int64)
    distances = Dict{Int64, Int64}()
    for vertex in g.vertices
        distances[vertex] = -1  # Use -1 to indicate unvisited
    end
    distances[v] = 0

    queue = Queue{Int64}()
    enqueue!(queue, v)

    while !isempty(queue)
        current = dequeue!(queue)
        current_distance = distances[current]
        for neighbor in get_neighbors(g, current)
            if distances[neighbor] == -1
                distances[neighbor] = current_distance + 1
                enqueue!(queue, neighbor)
            end
        end
    end
    
    return distances
end

function pairs_by_distance(g::QGraph)
    dist_pairs = Dict{Int, Array{Tuple{Int64, Int64}, 1}}()

    for v in g.vertices
        dist_from_v = distance(g, v)
        for (u, d) in dist_from_v
            if d > 0
                pair = v < u ? (v, u) : (u, v)
                if !haskey(dist_pairs, d)
                    dist_pairs[d] = Array{Tuple{Int64, Int64}, 1}()
                end
                if !(pair in dist_pairs[d])
                    push!(dist_pairs[d], pair)
                end
            end
        end
    end
    
    return dist_pairs
end