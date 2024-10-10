include("graph.jl")
include("distribution.jl")
include("topology.jl")
include("quantumrules.jl")

using ProgressBars

function less(a, b)
    return a < b
end

function less_or_equal(a, b)
    return a ≤ b
end

function no_comp(a, b)
    return true
end

function choose_edge(g::QGraph, v_start::Int64, v_end::Int64; ignore_nodes=[], compare=less, σ = 0.0)
    dist = distance(g, v_end)
    solutions = []
    λs = []
    Ns = []
    distances = []
    nodes = []

    if are_neighbors(g, v_start, v_end)
        push!(solutions, [])
        push!(λs, get_weights(g, v_start, v_end))
        push!(Ns, 0)
        push!(distances, 0)
        push!(nodes, v_end)
    end

    # First, operate on distance d = 0
    # Iterate for each neighbor
    for v_neighbor in filter(v -> compare(dist[v], dist[v_start]) && v ∉ ignore_nodes, get_neighbors(g, v_start))
        K_τ = get_weights(g, v_start, v_neighbor)

        tot_sol = []
        λ = K_τ

        N = 0
        # Check for distillations between v_start and v_neighbor
        if length(K_τ) > 1
            sort!(K_τ)
            λ = K_τ[1]
            
            while λ > 0.5 && N + 2 ≤ length(K_τ)
                λ = λdistill(λ, K_τ[N+2])
                push!(tot_sol, [v_start, v_neighbor])
                N += 1
            end
        end

        # If λ is not maximally entangled, look for swap solutions
        if λ > 0.5
            # Find all common neighbors between v and v_neighbor 
            swap_routes = filter(v -> v ≠ v_neighbor && v ∉ ignore_nodes && v_neighbor ∈ get_neighbors(g, v), get_neighbors(g, v_start))

            # Compute Schmidt values after swapping
            λ_swap_routes = Dict()
            for v in swap_routes
                λ_swap_routes[v] = λswap(get_weights(g, v_start, v), get_weights(g, v, v_neighbor))
            end

            # Sort routes by their predicted Schmidt value after swapping
            sort!(swap_routes, by = v -> λ_swap_routes[v])

            # Swap and distill until necessary
            i = 1
            while λ > 0.5 && i ≤ length(swap_routes)
                v = swap_routes[i]
                λ = λdistill(λ, λ_swap_routes[v])
                push!(tot_sol, [v_start, v, v_neighbor], [v_start, v_neighbor])
                N += 2
                i += 1
            end
        end
        if v_neighbor == v_end && λ == 0.5
            return tot_sol, v_end, 0.5
        end
        if !isempty(tot_sol)
            push!(solutions, tot_sol)
            push!(λs, λ)
            push!(Ns, N)
            push!(distances, dist[v_neighbor])
            push!(nodes, v_neighbor)
        end
    end

    # Now, check nodes at distance d = 1
    # Find nodes at distance d = 1
    v_distance_1 = Dict{Int64, Vector{Int64}}()
    vs = []
    for v_neighbor in filter(v_neigh -> v_neigh ∉ ignore_nodes, get_neighbors(g, v_start))
        v_distance_1[v_neighbor] = filter(v -> v ≠ v_start && v ∉ ignore_nodes &&
                                            compare(dist[v], dist[v_start]) && 
                                            !are_neighbors(g, v_start, v), get_neighbors(g, v_neighbor))
        append!(vs, v_distance_1[v_neighbor])
    end
    unique!(vs)

    # Iterate for each node at distance d = 1
    for v in vs
        tot_sol = []

        # Find all common neighbors of v_start and v
        swap_routes = [node for node in keys(v_distance_1) if v in v_distance_1[node]]

        # Compute the predicted Schmidt value after swapping
        λ_swap_routes = Dict()
        for v_neighbor in swap_routes
            λ_swap_routes[v_neighbor] = λswap(get_weights(g, v_start, v_neighbor), get_weights(g, v_neighbor, v))
        end
        sort!(swap_routes, by = node -> λ_swap_routes[node])

        # Swap the pair with predicted minimum λ
        λ = λ_swap_routes[swap_routes[1]]
        push!(tot_sol, [v_start, swap_routes[1], v])

        i = 2
        N = 1
        # Distill new swap links until maximum entanglement or until possible
        while λ > 0.5 && i ≤ length(swap_routes)
            v_neighbor = swap_routes[i]
            λ = λdistill(λ, λ_swap_routes[v_neighbor])
            push!(tot_sol, [v_start, v_neighbor, v], [v_start, v])
            N += 2
            i += 1
        end
        if v == v_end && λ == 0.5
            return tot_sol, v_end, 0.5
        end
        if !isempty(tot_sol)
            push!(solutions, tot_sol)
            push!(λs, λ)
            push!(Ns, N)
            push!(distances, dist[v])
            push!(nodes, v)
        end
    end

    if isempty(solutions)
        return nothing, nothing, nothing
    end
    λmin = minimum(λs)

    valid_sols = []
    for (i, λ) in enumerate(λs)
        if λ <= λmin + σ
            push!(valid_sols, i)
        end
    end

    selection = rand(valid_sols)

    return solutions[selection], nodes[selection], λs[selection]
end

function filter_path(paths, g::QGraph, v_start::Int, v_end::Int)
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

function find_path(g::QGraph, v_start::Int64, v_end::Int64; samples=64, σ=0.0, λmean=0.5, progressbar=true)
    start_end_dist = distance(g, v_end)[v_start]
    actual_σs = repeat(collect(range(0.0, σ, samples))[1:3:end], inner=3)[1:samples]

    paths = fill([], samples)
    λs = fill(0.0, samples)
    destroyed_links = fill([], samples)

    Threads.@threads for i in (progressbar ? ProgressBar(1:samples) : 1:samples)
    # for i in (progressbar ? ProgressBar(1:samples) : 1:samples)
        g_new = deepcopy(g)
        new_path = []

        compare = no_comp
        if i % 3 == 1
            compare = less
        elseif i % 3 == 2
            compare = less_or_equal
        end

        if are_neighbors(g, v_start, v_end)
            if get_weights(g, v_start, v_end) > 0.5
                new_path, node, old_λ = choose_edge(g_new, v_start, v_end, σ=actual_σs[i])
                if node == v_end
                    paths[i] = new_path
                    λs[i] = old_λ
                else
                    old_λ = get_weights(g, v_start, v_end)
                end
                if old_λ > 0.5
                    max_dist_search = 10
                    sample_dist = 9
                    σs_dist = repeat(collect(range(0.0, σ, sample_dist))[1:3:end], inner=3)[1:sample_dist]

                    n1, n2 = v_start, v_end
                    for s in 1:sample_dist
                        # println("Looking for new path between $n1 and $n2")
                        g_new_new = deepcopy(g)
                        rem_edge!(g_new_new, v_start, v_end, get_weights(g, v_start, v_end))

                        compare_dist = no_comp
                        if s % 3 == 1
                            compare_dist = less
                        elseif s % 3 == 2
                            compare_dist = less_or_equal
                        end

                        iter_dist = 0
                        new_n1 = n1
                        alternative_path = []
                        while !are_neighbors(g_new_new, n1, n2) && iter_dist < max_dist_search
                            if iter_dist > div(max_dist_search, 2) && compare_dist == no_comp
                                compare_dist = less
                            end
                            new_solution, next_node_dist = choose_edge(g_new_new, new_n1, n2, ignore_nodes=unique(vcat(paths[i]...)), compare=compare_dist, σ=σs_dist[s])
                            if !isnothing(new_solution)
                                if new_n1 ≠ n1
                                    push!(new_solution, [n1, new_n1, next_node_dist])
                                end
                                final_new_sol = []
                                for (k, op) in enumerate(new_solution)
                                    if length(op) == 3
                                        push!(final_new_sol, op)
                                        entanglement_swapping!(g_new_new, op[1], op[2], op[3])
                                        if length(get_weights(g_new_new, op[1], op[3])) > 1 && ((k < length(new_solution) && sort(new_solution[k+1]) != sort([op[1], op[3]])) || k == length(new_solution))
                                            push!(final_new_sol, [op[1], op[3]])
                                            entanglement_distillation!(g_new_new, op[1], op[3])
                                        end
                                    elseif length(op) == 2
                                        push!(final_new_sol, op)
                                        entanglement_distillation!(g_new_new, op[1], op[2])
                                    end
                                end
                                # println("Found solution: ", final_new_sol)
                                append!(alternative_path, final_new_sol)
                                new_solution = final_new_sol
                                new_n1 = next_node_dist
                            end
                            iter_dist += isnothing(new_solution) ? max_dist_search : 1
                        end

                        if iter_dist < max_dist_search && !isempty(alternative_path) && 
                                    λdistill(get_weights(g_new_new, n1, n2), old_λ) < old_λ
                            append!(paths[i], alternative_path)
                            push!(paths[i], [v_start, v_end])
                            old_λ = λdistill(get_weights(g_new_new, n1, n2), old_λ)
                            if old_λ == 0.5
                                break
                            end
                        end
                    end
                    if isempty(paths[i])
                        λs[i] = get_weights(g, v_start, v_end)
                    else
                        g_new = apply(g, paths[i])
                        λs[i] = get_weights(g_new, v_start, v_end)
                        destroyed_links[i] = [[v1, v2] for v1 in unique(collect(Iterators.flatten(paths[i]))) 
                                                    for v2 in unique(collect(Iterators.flatten(paths[i])))
                                                    if v1 < v2 && are_neighbors(g, v1, v2) && (!are_neighbors(g_new, v1, v2) || 
                                                        ((v1, v2) == tuple(sort([v_start, v_end])...) && !isempty(paths[i][1])))]
                    end
                else
                    destroyed_links[i] = [[v1, v2] for v1 in unique(collect(Iterators.flatten(paths[i]))) 
                                                    for v2 in unique(collect(Iterators.flatten(paths[i])))
                                                    if v1 < v2 && are_neighbors(g, v1, v2) && (!are_neighbors(g_new, v1, v2) || 
                                                        ((v1, v2) == tuple(sort([v_start, v_end])...) && !isempty(paths[i][1])))]
                end
            else
                λs[i] = 0.5
            end
        else
            iter = 0
            explored_nodes = []
            explored_λs = []
            all_paths = []
            new_v_start = v_start
            while !are_neighbors(g_new, v_start, v_end) && iter < 4*start_end_dist
                solution, next_node, next_λ = choose_edge(g_new, new_v_start, v_end, ignore_nodes=[v_start], compare=compare, σ=actual_σs[i])
                if !isnothing(solution)
                    push!(explored_nodes, next_node)
                    if new_v_start ≠ v_start
                        push!(solution, [v_start, new_v_start, next_node])
                    end
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

                    push!(all_paths, solution)
                    push!(explored_λs, next_λ)
                    append!(new_path, solution)
                end
                iter += isnothing(solution) ? 4*start_end_dist : 1
                new_v_start = next_node
            end

            if iter < 4*start_end_dist
                if get_weights(g_new, v_start, v_end) > 0.5
                    max_dist_search = 10
                    # println("Previous path: ", new_path)
                    explored_nodes = [v_start, explored_nodes...]
                    for (idx1, idx2) in [(i-1, i) for i in 2:length(explored_nodes) if explored_λs[i-1] > 0.5]
                        n1, n2 = explored_nodes[idx1], explored_nodes[idx2]
                        # println("Looking for new path between $n1 and $n2")
                        g_new_new = deepcopy(g)
                        for i in 1:idx1
                            g_new_new = apply(g_new_new, all_paths[i])
                        end
                        

                        iter_dist = 0
                        new_n1 = n1
                        alternative_path = []
                        while !are_neighbors(g_new_new, n1, n2) && iter_dist < max_dist_search
                            new_solution, next_node_dist = choose_edge(g_new_new, new_n1, n2, ignore_nodes=unique(vcat(new_path...)), compare=(iter_dist < div(max_dist_search, 2) ? no_comp : less), σ=actual_σs[i])
                            if !isnothing(new_solution)
                                if new_n1 ≠ n1
                                    push!(new_solution, [n1, new_n1, next_node_dist])
                                end
                                final_new_sol = []
                                for (k, op) in enumerate(new_solution)
                                    if length(op) == 3
                                        push!(final_new_sol, op)
                                        entanglement_swapping!(g_new_new, op[1], op[2], op[3])
                                        if length(get_weights(g_new_new, op[1], op[3])) > 1 && ((k < length(new_solution) && sort(new_solution[k+1]) != sort([op[1], op[3]])) || k == length(new_solution))
                                            push!(final_new_sol, [op[1], op[3]])
                                            entanglement_distillation!(g_new_new, op[1], op[3])
                                        end
                                    elseif length(op) == 2
                                        push!(final_new_sol, op)
                                        entanglement_distillation!(g_new_new, op[1], op[2])
                                    end
                                end
                                # println("Found solution: ", final_new_sol)
                                append!(alternative_path, final_new_sol)
                                new_solution = final_new_sol
                                new_n1 = next_node_dist
                            end
                            iter_dist += isnothing(new_solution) ? max_dist_search : 1
                        end

                        if iter_dist < max_dist_search && !isempty(alternative_path)
                            # println("Found path: ", alternative_path)
                            final_new_path = []
                            for i in 2:idx1
                                append!(final_new_path, all_paths[i-1])
                            end
                            idx_add_later = n1 == v_start ? length(all_paths[idx1])+1 : findfirst(op -> op == [v_start, n1, n2], all_paths[idx1])
                            append!(final_new_path, all_paths[idx1][1:idx_add_later-1])
                            append!(final_new_path, alternative_path)

                            push!(final_new_path, [n1, n2])
                            append!(final_new_path, all_paths[idx1][idx_add_later:end])
                            for i in idx2:length(all_paths)
                                append!(final_new_path, all_paths[i])
                            end
                            new_path = final_new_path
                            # println("Final path: $(new_path)")
                        end
                    end

                    # println()
                    # println("Final path: $(new_path)")
                    # println()
                    # println()
                    # println()
                    # println()
                    # println()
                    # println()
                    # println()
                end
                
                g_new = apply(g, new_path)
                # new_path, g_new = filter_path(new_path, g, v_start, v_end)
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
        if λs[i] < λmin || (isapprox(λs[i], λmin, atol=1e-3) && length(paths[i]) < N)
            path = deepcopy(paths[i])
            λmin = λs[i]
            N = length(paths[i])
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
                N = length(paths[i]) + length(paths[j])
            end
        end
    end

    if isnothing(path)
        return nothing
    end

    return path, λmin, N
end