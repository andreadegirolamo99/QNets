include("src/percolation.jl")


function run_mean(n, i, λmean, topology)
    g = eval(Symbol(topology))(n, λ=λmean)

    pairs = pairs_by_distance(g)

    println("Mean of links: $λmean")
    flush(stdout)
    mean_ent = Dict()
    mean_N = Dict()

    for distance in sort(collect(keys(pairs)))[2:end]
        ent = fill(0.0, length(pairs[distance]))
        Ns = fill(0.0, length(pairs[distance]))
        println("Computing entanglement for nodes at distance $(distance-1), λ = $λmean")
        flush(stdout)

        Threads.@threads for j in 1:length(pairs[distance])
        # for j in 1:length(pairs[distance])
            v1, v2 = pairs[distance][j]
            sol = find_path(g, v1, v2, samples=600, σ=0.5, progressbar=false)
            if !isnothing(sol)
                path, λ, N = sol
                ent[j] = 2*(1 - λ)
                Ns[j] = N
            end
        end
        println()
        flush(stdout)

        mean_ent[distance] = mean(ent)
        mean_N[distance] = mean(Ns)
    end

    file = open("out_files/$(topology)_$n/mean_/$(i)_$(λmean).out", "w")
    for dist in sort(collect(keys(mean_ent)))
        write(file, "$(dist-1) $(mean_ent[dist]) $(mean_N[dist]) \n")
    end

    flush(file)
    close(file)
end

function run_dev(n, i, λmean, σ, topology, sample)
    g = eval(Symbol(topology))(n, σ=σ, λmean=λmean)

    pairs = pairs_by_distance(g)

    @assert isapprox(mean(get_weights(g, pair...) for pair in pairs[1]), λmean)

    println("Mean of links: $λmean")
    flush(stdout)
    mean_ent = Dict()
    mean_N = Dict()
    
    for distance in sort(collect(keys(pairs)))[2:end]
        ent = fill(0.0, length(pairs[distance]))
        Ns = fill(0.0, length(pairs[distance]))
        println("Computing entanglement for nodes at distance $(distance-1), λ = $λmean")
        flush(stdout)

        Threads.@threads for j in 1:length(pairs[distance])
            v1, v2 = pairs[distance][j]
            sol = find_path(g, v1, v2, samples=600, σ=0.5, progressbar=false)
            if !isnothing(sol)
                path, λ, N = sol
                ent[j] = 2*(1 - λ)
                Ns[j] = N
            end
        end
        println()
        flush(stdout)

        mean_ent[distance] = mean(ent)
        mean_N[distance] = mean(Ns)
    end

    file = open("out_files/$(replace(topology, r"noisy_" => ""))_$n/std_dev=$(σ)/$(sample)/$(i)_$(λmean).out", "w")
    for dist in sort(collect(keys(mean_ent)))
        write(file, "$(dist-1) $(mean_ent[dist]) $(mean_N[dist]) \n")
    end
    flush(file)
    close(file)
end

function run(topology, N, i, λmean, std_dev, sample)
    if std_dev == 0
        run_mean(N, i, λmean, topology)
    else
        run_dev(N, i, λmean, std_dev, topology, sample)
    end
end

topology = ARGS[1]
N = parse(Int64, ARGS[2])
i = parse(Int64, ARGS[3])
λmean = parse(Float64, ARGS[4])
std_dev = parse(Float64, ARGS[5])
sample = parse(Int64, ARGS[6])

run(topology, N, i, λmean, std_dev, sample)

# run_dev(6, 1, 0.01)