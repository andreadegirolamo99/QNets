include("src/percolation.jl")


function percolate(g)
    pairs = pairs_by_distance(g)

    # Cache sorted distances
    sorted_distances = sort(collect(keys(pairs)))[2:end]
    
    mean_ent = fill(0.0, length(sorted_distances))
    mean_N = fill(0.0, length(sorted_distances))

    # Parallelize across all distances using spawn
    results = Vector{Task}(undef, length(sorted_distances))

    for idx in eachindex(sorted_distances)
        results[idx] = Threads.@spawn begin
            distance = sorted_distances[idx]
            ent = Vector{Float64}(undef, length(pairs[distance]))
            Ns = Vector{Float64}(undef, length(pairs[distance]))

            # Process pairs in parallel
            Threads.@threads for j in 1:length(pairs[distance])
                v1, v2 = pairs[distance][j]
                sol = find_path(g, v1, v2, samples=600, σ=0.5, progressbar=false)
                if !isnothing(sol)
                    path, λ, N = sol
                    ent[j] = 2 * (1 - λ)
                    Ns[j] = N
                else
                    ent[j] = 0.0  # Handle case where no solution is found
                    Ns[j] = 0.0
                end
            end

            println("Finished distance $(distance - 1)")
            flush(stdout)

            return mean(ent), mean(Ns)
        end
    end

    # Wait for all tasks to complete and collect the results
    for idx in 1:length(sorted_distances)
        mean_ent[idx], mean_N[idx] = fetch(results[idx])
    end

    return sorted_distances, mean_ent, mean_N
end

function run(topology, N, i, λmean, std_dev, sample)
    println("Mean of links: $λmean")
    println()
    flush(stdout)

    if std_dev == 0
        g = eval(Symbol(topology))(N, λ=λmean)
        sorted_distances, mean_ent, mean_N = percolate(g)

        file = open("out_files/$(topology)_$N/mean/$(i)_$(λmean).out", "w")
        
        for dist_idx in eachindex(mean_ent)
            write(file, "$(sorted_distances[dist_idx]-1) $(mean_ent[dist_idx]) $(mean_N[dist_idx]) \n")
        end

        flush(file)
        close(file)
    else
        g = eval(Symbol(topology))(N, σ=σ, λmean=λmean)
        @assert isapprox(mean(get_weights(g, pair...) for pair in pairs[1]), λmean)
        sorted_distances, mean_ent, mean_N = percolate(g)

        file = open("out_files/$(replace(topology, r"noisy_" => ""))_$N/std_dev=$(σ)/$(sample)/$(i)_$(λmean).out", "w")
        for dist_idx in eachindex(mean_ent)
            write(file, "$(sorted_distances[dist_idx]-1) $(mean_ent[dist_idx]) $(mean_N[dist_idx]) \n")
        end
        flush(file)
        close(file)
    end
end


topology = ARGS[1]
N = parse(Int64, ARGS[2])
i = parse(Int64, ARGS[3])
λmean = round(parse(Float64, ARGS[4]), digits=3)
std_dev = round(parse(Float64, ARGS[5]), digits=2)
sample = parse(Int64, ARGS[6])

run(topology, N, i, λmean, std_dev, sample)

# run_dev(6, 1, 0.01)