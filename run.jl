include("src/QSN.jl")

using .QSN
using Distributions
using SparseArrays
using Combinatorics
using LinearAlgebra
using ITensors
using DataStructures
using Plots
using GraphPlot
using GraphRecipes
using ProgressBars

P(λ, δ) = 2*(2*(1 - λ))^(1/δ - 1)/δ

sample_λ(δ) = round(0.5*(rand(Beta(1, 1/δ)) + 1), digits=2)

mutable struct Graph
    vertices::Array{Int64}
    edges::Dict{Tuple{Int64, Int64}, Union{Float64, Array{Float64, 1}}}

    Graph() = new([], Dict())

    Graph(n::Int64) = new([i for i in 1:n], Dict())
end

add_vertex!(g::Graph) = push!(g.vertices, length(g.vertices)+1)
add_vertices!(g::Graph, n::Int64) = append!(g.vertices, collect(length(g.vertices)+1:length(g.vertices)+1+n))

vertices(g::Graph) = g.vertices
n_vertices(g::Graph) = length(g.vertices)

function add_edge!(g::Graph, v1::Int64, v2::Int64, w::Float64)
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

get_weights(g::Graph, v1::Int64, v2::Int64) = are_neighbors(g, v1, v2) ? (v1 < v2 ? g.edges[(v1, v2)] : g.edges[(v2, v1)]) : []
    
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

edges(g::Graph, v::Int64) = filter(edge -> v in edge, keys(g.edges))

function get_neighbors(g::Graph, v::Int64)
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

are_neighbors(g::Graph, v1::Int64, v2::Int64) = v2 in get_neighbors(g, v1)

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

function adjMat(g::Graph)
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

function square_lattice(n::Int64; λ=0.5)
    g = Graph(n^2)

    for j in 1:n
        for i in 1:n-1
            add_edge!(g, i+n*(j-1), i+n*(j-1)+1, λ)
        end
    end
    for j in 1:n-1
        for i in 1:n
            add_edge!(g, i+n*(j-1), i+n*j, λ)
        end
    end

    return g
end
function diagonal_square_lattice(n::Int64; λ=0.5)
    g = Graph(n^2)

    for j in 1:n-1
        for i in 1:n-1
            add_edge!(g, i+n*(j-1), i+n*(j-1)+1, λ)
            add_edge!(g, i+n*(j-1), i+n*j+1, λ)
        end
    end
    for i in 1:n-1
        add_edge!(g, i+n*(n-1), i+n*(n-1)+1, λ)
    end

    for j in 1:n-1
        for i in 1:n-1
            add_edge!(g, i+n*(j-1), i+n*j, λ)
            add_edge!(g, i+n*(j-1)+1, i+n*j, λ)
        end
        add_edge!(g, n*j, n*(j+1), λ)
    end

    return g
end

function noisy_square_lattice(n::Int64; δ=1e-6)
    g = Graph(n^2)

    for j in 1:n
        for i in 1:n-1
            add_edge!(g, i+n*(j-1), i+n*(j-1)+1, sample_λ(δ))
        end
    end
    for j in 1:n-1
        for i in 1:n
            add_edge!(g, i+n*(j-1), i+n*j, sample_λ(δ))
        end
    end

    return g
end

function noisy_diagonal_square_lattice(n::Int64; δ=1e-6)
    g = Graph(n^2)

    for j in 1:n-1
        for i in 1:n-1
            add_edge!(g, i+n*(j-1), i+n*(j-1)+1, sample_λ(δ))
            add_edge!(g, i+n*(j-1), i+n*j+1, sample_λ(δ))
        end
    end
    for i in 1:n-1
        add_edge!(g, i+n*(n-1), i+n*(n-1)+1, sample_λ(δ))
    end

    for j in 1:n-1
        for i in 1:n-1
            add_edge!(g, i+n*(j-1), i+n*j, sample_λ(δ))
            add_edge!(g, i+n*(j-1)+1, i+n*j, sample_λ(δ))
        end
        add_edge!(g, n*j, n*(j+1), sample_λ(δ))
    end

    return g
end

λswap(α, β) = 0.5*(1 + √(1 - 16*α*(1-α)*β*(1-β)))

function entanglement_swapping!(g, v1, v2, v3)
    α = get_weights(g, v1, v2)
    β = get_weights(g, v2, v3)

    λ = λswap(α, β)

    rem_edge!(g, min(v1, v2), max(v1, v2), α)
    rem_edge!(g, min(v2, v3), max(v2, v3), β)
    add_edge!(g, min(v1, v3), max(v1, v3), λ)
end

λdistill(α, β) = max(0.5, α*β)

function entanglement_distillation!(g, v1, v2)
    α, β = get_weights(g, v1, v2)[1:2]

    λ = λdistill(α, β)

    rem_edge!(g, min(v1, v2), max(v1, v2), α)
    rem_edge!(g, min(v1, v2), max(v1, v2), β)
    add_edge!(g, min(v1, v2), max(v1, v2), λ)
end

function connect_ends(g::Graph)
    n = Int64(√(n_vertices(g)))

    for i in 1:n-1
        entanglement_swapping!(g, n*(i - 1) + i, n*(i - 1) + i + 1, (n + 1)*i + 1)
        entanglement_swapping!(g, n*(i - 1) + i, (n + 1)*i, (n + 1)*i + 1)
        entanglement_distillation!(g, n*(i - 1) + i, (n + 1)*i + 1)
    end

    for i in 1:n-2
        entanglement_swapping!(g, 1, (n + 1)*i + 1, (n + 1)*(i + 1) + 1)
    end

    return get_weights(g, 1, n_vertices(g))
end

function distance(g::Graph, v::Int64)
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

function distillation_sols(close_neighbors::Dict{Int64, Vector{Int64}}; exclude=[])
    node_to_keys = Dict{Int64, Vector{Int64}}()

    for (key, nodes) in close_neighbors
        for node in nodes
            if !haskey(node_to_keys, node)
                node_to_keys[node] = Int[]
            end
            push!(node_to_keys[node], key)
        end
    end

    return sort.(collect.(Set(Set.([[[key1, node], [key2, node]] for (node, keys) in node_to_keys 
                                        if length(keys) > 1 for key1 in keys for key2 in keys 
                                            if key1 ≠ key2 && key1 ∉ exclude && key2 ∉ exclude && node ∉ exclude]))))
end

function choose_edge(g::Graph, v_start::Int64, v_end::Int64; σ = 0.0, λmean = 0.5)
    dist = distance(g, v_end)
    solutions = []
    λs = []
    Ns = []

    if are_neighbors(g, v_start, v_end)
        push!(solutions, [])
        push!(λs, get_weights(g, v_start, v_end))
        push!(Ns, 0)
    end

    close_neighbors = Dict{Int64, Vector{Int64}}()
    for v in filter(v -> dist[v] ≤ dist[v_start], get_neighbors(g, v_start))
        close_neighbors[v] = filter(n -> n ≠ v_start && dist[n] ≤ dist[v_start], get_neighbors(g, v))
    end

    for v1 in keys(close_neighbors)
        for v2 in close_neighbors[v1]
            α = get_weights(g, v_start, v1)
            β = get_weights(g, v1, v2)

            push!(solutions, [[v_start, v1, v2]])
            push!(λs, λswap(α, β))
            push!(Ns, 2)
        end
    end

    if λmean ≠ 0.5
        for distill in distillation_sols(close_neighbors)
            v1 = distill[1][1]
            v2 = distill[2][1]
            v3 = distill[2][2]
    
            α1 = get_weights(g, v_start, v1)
            β1 = get_weights(g, v1, v3)
            α2 = get_weights(g, v_start, v2)
            β2 = get_weights(g, v2, v3)
    
            push!(solutions, [[v_start, v1, v3], [v_start, v2, v3], [v_start, v3]])
            push!(λs, λdistill(λswap(α1, β1), λswap(α2, β2)))
            push!(Ns, 4)
        end
    
        for v1 in get_neighbors(g, v_start)
            close_neighbors_of_neighbor = Dict{Int64, Vector{Int64}}()
            for v2 in filter(v2 -> v2 ≠ v_start &&  dist[v2] ≤ dist[v_start], get_neighbors(g, v1))
                close_neighbors_of_neighbor[v2] = filter(v3 -> v3 ≠ v1 && dist[v3] ≤ dist[v_start], get_neighbors(g, v2))
            end
    
            for distill_neighbor in distillation_sols(close_neighbors_of_neighbor, exclude=[v_start])
                v2 = distill_neighbor[1][1]
                v3 = distill_neighbor[2][1]
                v4 = distill_neighbor[2][2]
    
                α = get_weights(g, v_start, v1)
                α1 = get_weights(g, v1, v2)
                β1 = get_weights(g, v2, v4)
                α2 = get_weights(g, v1, v3)
                β2 = get_weights(g, v3, v4)
    
                push!(solutions, [[v1, v2, v4], [v1, v3, v4], [v1, v4], [v_start, v1, v4]])
                push!(λs, λswap(α, λdistill(λswap(α1, β1), λswap(α2, β2))))
                push!(Ns, 4)
            end
        end

        for v1 in filter(v -> dist[v] ≤ dist[v_start], get_neighbors(g, v_start))
            for v_common in intersect(close_neighbors[v1], get_neighbors(g, v_start))
                α1 = get_weights(g, v_start, v1)
                α2 = get_weights(g, v1, v_common)
                α3 = get_weights(g, v_start, v_common)

                push!(solutions, [[v_start, v1, v_common], [v_start, v_common]])
                push!(λs, λdistill(λswap(α1, α2), α3))
                push!(Ns, 2)
            end
        end
    end

    if isempty(solutions)
        return nothing
    end
    λmin = minimum(λs)

    valid_sols = []
    for (i, λ) in enumerate(λs)
        if λ <= λmin + σ
            push!(valid_sols, i)
        end
    end

    return solutions[rand(valid_sols)]
end

function filter_path(paths, g::Graph, v_start::Int, v_end::Int)
    connections = [[v] for v in get_neighbors(g, v_start)]
    pending = []
    ignore = false

    for i in eachindex(paths)
        op = paths[i]
        if length(op) == 3
            if ignore
                ignore = false
            else
                if op[1] == v_start
                    push!(connections[findfirst(seq -> seq[end] == op[2], connections)], op[3])
                else
                    push!(pending, op)
                end
            end
        else
            if !isempty(pending)
                if length(pending) > 1
                    push!(connections[findfirst(seq -> seq[end] == pending[1][1], connections)], pending[1][2], pending[2][2], pending[2][3])
                    pending = []
                    ignore = true
                end
            else
                idx = op[1] == v_start ? 2 : 1
                paths_to_merge = findall(seq -> seq[end] == op[idx], connections)
                connections[paths_to_merge[1]] = filter(v -> v != op[idx], unique([connections[paths_to_merge[1]]..., connections[paths_to_merge[2]]...]))
                deleteat!(connections, paths_to_merge[2])
                push!(connections[paths_to_merge[1]], op[idx])
            end
        end
    end

    sequence = [v_start, connections[findfirst(seq -> seq[end] == v_end, connections)]...]
    g_new = deepcopy(g)

    final_path = []
    for i in eachindex(paths)
        path = paths[i]
        if all(node in sequence for node in path)
            if length(path) == 3 && are_neighbors(g_new, path[1], path[2]) && are_neighbors(g_new, path[2], path[3])
                push!(final_path, path)
                entanglement_swapping!(g_new, path...)
                if length(get_weights(g_new, path[1], path[3])) > 1 && ((i+1 <= length(paths) && sort(paths[i+1]) != sort([path[1], path[3]])) || i == length(paths))
                    push!(final_path, [path[1], path[3]])
                    entanglement_distillation!(g_new, path[1], path[3])
                end
            else
                if are_neighbors(g_new, path[1], path[2]) && length(get_weights(g_new, path[1], path[2])) > 1
                    push!(final_path, path)
                    entanglement_distillation!(g_new, path...)
                    if sort([path[1], path[2]]) == sort([v_start, v_end]) && get_weights(g, v_start, v_end) == 0.5
                        return final_path
                    end
                end
            end
        end
    end

    return final_path, g_new
end

function find_path(g::Graph, v_start::Int64, v_end::Int64; samples=64, σ=0.0, λmean=0.5, progressbar=true)

    start_end_dist = distance(g, v_end)[v_start]
    actual_σs = collect(range(0.0, σ, samples))

    unfiltered_paths = fill([], samples)
    paths = fill([], samples)
    λs = fill(0.0, samples)
    destroyed_links = fill([], samples)

    Threads.@threads for i in (progressbar ? ProgressBar(1:samples) : 1:samples)
        g_new = deepcopy(g)
        new_path = []

        if are_neighbors(g, v_start, v_end)
            if get_weights(g, v_start, v_end) == 0.5
                λs[i] = 0.5
            else
                iter = 0
                prev_λ = get_weights(g, v_start, v_end)
                new_connection = false
                while !new_connection && iter < 2*length(vertices(g))
                    solution = choose_edge(g_new, v_start, v_end, σ=actual_σs[i], λmean=λmean)
                    if !isnothing(solution)
                        final_sol = []
                        for (k, op) in enumerate(solution)
                            if length(op) == 3
                                push!(final_sol, op)
                                entanglement_swapping!(g_new, op[1], op[2], op[3])
                                if length(get_weights(g_new, op[1], op[3])) > 1 && ((k < length(solution) && sort(solution[k+1]) != sort([op[1], op[3]])) || k == length(solution))
                                    push!(final_sol, [op[1], op[3]])
                                    entanglement_distillation!(g_new, op[1], op[3])
                                end
                            elseif length(op) == 2
                                push!(final_sol, op)
                                entanglement_distillation!(g_new, op[1], op[2])
                            end
                        end

                        solution = final_sol

                        distills = [[v1, v2] for v1 in unique(collect(Iterators.flatten(solution))) 
                                            for v2 in unique(collect(Iterators.flatten(solution)))
                                            if v1 < v2 && are_neighbors(g_new, v1, v2) && length(get_weights(g_new, v1, v2)) == 2]
                                                
                        if !isempty(distills)
                            distills = first(distills)
                            entanglement_distillation!(g_new, distills[1], distills[2])
                            push!(solution, distills)
                        end
                        append!(new_path, solution)
                    end
                    new_connection = isnothing(solution) || (are_neighbors(g_new, v_start, v_end) && 
                                                            get_weights(g_new, v_start, v_end) < prev_λ &&
                                                            !isapprox(get_weights(g_new, v_start, v_end), prev_λ))
                    iter += isnothing(solution) ? 2*length(vertices(g)) : 1
                end

                if iter < 2*length(vertices(g))
                    unfiltered_paths[i] = new_path
                    new_path, g_new = filter_path(new_path, g, v_start, v_end)
                    paths[i] = new_path
                    λs[i] = get_weights(g_new, v_start, v_end)
                    destroyed_links[i] = [[v1, v2] for v1 in unique(collect(Iterators.flatten(new_path))) 
                                                for v2 in unique(collect(Iterators.flatten(new_path)))
                                                if v1 < v2 && are_neighbors(g, v1, v2) && (!are_neighbors(g_new, v1, v2) || 
                                                    ((v1, v2) == tuple(sort([v_start, v_end])...) && !isempty(paths[i][1])))]
                else
                    λs[i] = prev_λ
                end
            end
        else
            iter = 0
            while !are_neighbors(g_new, v_start, v_end) && iter < 4*start_end_dist
                solution = choose_edge(g_new, v_start, v_end, σ=actual_σs[i], λmean=λmean)
                if !isnothing(solution)
                    final_sol = []
                    for (k, op) in enumerate(solution)
                        if length(op) == 3
                            push!(final_sol, op)
                            entanglement_swapping!(g_new, op[1], op[2], op[3])
                            if length(get_weights(g_new, op[1], op[3])) > 1 && ((k < length(solution) && sort(solution[k+1]) != sort([op[1], op[3]])) || k == length(solution))
                                push!(final_sol, [op[1], op[3]])
                                entanglement_distillation!(g_new, op[1], op[3])
                            end
                        elseif length(op) == 2
                            push!(final_sol, op)
                            entanglement_distillation!(g_new, op[1], op[2])
                        end
                    end

                    solution = final_sol

                    distills = [[v1, v2] for v1 in unique(collect(Iterators.flatten(solution))) 
                                            for v2 in unique(collect(Iterators.flatten(solution)))
                                            if v1 < v2 && are_neighbors(g_new, v1, v2) && length(get_weights(g_new, v1, v2)) == 2]
        
                    if !isempty(distills)
                        distills = first(distills)
                        entanglement_distillation!(g_new, distills[1], distills[2])
                        push!(solution, distills)
                    end

                    append!(new_path, solution)
                end
                iter += isnothing(solution) ? 4*start_end_dist : 1
            end

            if iter < 4*start_end_dist
                unfiltered_paths[i] = new_path
                new_path, g_new = filter_path(new_path, g, v_start, v_end)
                paths[i] = new_path
                λs[i] = get_weights(g_new, v_start, v_end)
                destroyed_links[i] = [[v1, v2] for v1 in unique(collect(Iterators.flatten(new_path))) 
                                                for v2 in unique(collect(Iterators.flatten(new_path)))
                                                if v1 < v2 && are_neighbors(g, v1, v2) && (!are_neighbors(g_new, v1, v2) || 
                                                    ((v1, v2) == tuple(sort([v_start, v_end])...) && !isempty(paths[i][1])))]
            else
                paths[i] = ["NOT FOUND"]
            end
        end
    end

    λmin = 1.0
    path = nothing
    N = Inf

    for i in findall(path -> path ≠ ["NOT FOUND"], paths)
        # println("Path to connect $(v_start) and $(v_end) before and after filtering:")
        # println(unfiltered_paths[i])
        # println(paths[i])
        # println("Entanglement: $(2*(1-λs[i]))")
        # println("Destroyed links: $(length(destroyed_links[i]))")
        # println()
        # println()
        # println()
        if λs[i] < λmin || (isapprox(λs[i], λmin, atol=1e-3) && length(destroyed_links[i]) < N)
            path = deepcopy(paths[i])
            λmin = λs[i]
            N = length(destroyed_links[i])
        end
        
        for j in 1:i-1
            independent = !isempty(paths[i]) && !isempty(paths[j])
            for link in destroyed_links[i]
                if link ∈ destroyed_links[j]
                    independent = false
                    break
                end
            end

            if independent && λs[i] > 0.5 && λs[j] > 0.5 && λmin > λdistill(λs[i], λs[j])
                path = [paths[i]..., paths[j]..., [v_start, v_end]]
                λmin = λdistill(λs[i], λs[j])
                N = length(destroyed_links[i]) + length(destroyed_links[j])
            end
        end
    end

    if isnothing(path)
        return nothing
    end

    return path, λmin, N
end

function apply(g::Graph, path)
    g_new = deepcopy(g)

    for op in path
        if length(op) == 3
            entanglement_swapping!(g_new, op[1], op[2], op[3])
        elseif length(op) == 2
            entanglement_distillation!(g_new, op[1], op[2])
        end
    end

    return g_new
end

function pairs_by_distance(g::Graph)
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

function run_pure(n)

    if !isdir("out_files/diag_$(n)")
        mkdir("out_files/diag_$(n)")
    end
    
    λs = collect(range(0.5, stop=1.0, length=100))

    for (i, λmean) in enumerate(λs)
        g = diagonal_square_lattice(n, λ=λmean)

        pairs = pairs_by_distance(g)

        println("Mean of links: $λmean")
        flush(stdout)
        mean_ent = Dict()
        mean_N = Dict()

        for distance in sort(collect(keys(pairs)))
            ent = []
            Ns = []
            println("Computing entanglement for nodes at distance $distance, λ = $λmean")
            flush(stdout)

            for (v1, v2) in (pairs[distance])
                sol = find_path(g, v1, v2, samples=1000, σ=0.15, λmean=λmean, progressbar=false)
                if !isnothing(sol)
                    path, λ, N = sol
                    push!(ent, 2*(1-λ))
                    push!(Ns, N)
                end
            end
            println()
            flush(stdout)

            mean_ent[distance] = mean(ent)
            mean_N[distance] = mean(Ns)
        end
        file = open("out_files/diag_$n/res_$(n)_$(i).out", "w")
        for dist in sort(collect(keys(mean_ent)))
            write(file, "$dist $(mean_ent[dist]) $(mean_N[dist]) \n")
        end
        flush(file)
        close(file)
    end
end

function run_mixed(n_samples, n)

    if !isdir("out_files/$(n)")
        mkdir("out_files/$(n)")
    end
    s_num = isempty(readdir("out_files/$n")) ? 1 : (maximum([parse(Int64, sample) for sample in readdir("out_files/$n")]) + 1)
    for i in s_num:s_num+n_samples-1
        if !isdir("out_files/$n/$i")
            mkdir("out_files/$n/$i")
        end
    end

    for j in 1:n_samples
    
        δs = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]#exp10.(range(-3, stop=0, length=12))
    
        for (i, δ) in enumerate(δs)
            g = noisy_square_lattice(n, δ=δ)
    
            pairs = pairs_by_distance(g)

            λmean = mean(get_weights(g, pair...) for pair in pairs[1])
            println("δ = $δ")
            println("Mean of links: $λmean")
            flush(stdout)
            mean_ent = Dict()
            mean_N = Dict()
    
            for distance in sort(collect(keys(pairs)))
                ent = []
                Ns = []
                println("Computing entanglement for nodes at distance $distance, δ = $δ")
                flush(stdout)
    
                for (v1, v2) in (pairs[distance])
                    sol = find_path(g, v1, v2, samples=1000, σ=0.15, λmean=λmean, progressbar=false)
                    if !isnothing(sol)
                        path, λ, N = sol
                        push!(ent, 2*(1-λ))
                        push!(Ns, N)
                    end
                end
                println()
                flush(stdout)
    
                mean_ent[distance] = mean(ent)
                mean_N[distance] = mean(Ns)
            end
            file = open("out_files/$n/$(s_num+j-1)/res_$(n)_$(i).out", "w")
            for dist in sort(collect(keys(mean_ent)))
                write(file, "$dist $(mean_ent[dist]) $(mean_N[dist]) \n")
            end
            flush(file)
            close(file)
        end
    end
end

# run_noise(16, 10)
run_pure(6)