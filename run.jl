include("src/percolation.jl")


function run_mean(n)

    if !isdir("out_files/diag_$(n)")
        mkdir("out_files/diag_$(n)")
    end
    
    λs = [0.6, 0.66, 0.7, 0.73]

    for (i, λmean) in enumerate(λs)
        g = diagonal_square_lattice(n, λ=λmean)

        pairs = pairs_by_distance(g)

        println("Mean of links: $λmean")
        flush(stdout)
        mean_ent = Dict()
        mean_N = Dict()

        for distance in sort(collect(keys(pairs)))[2:end]
            ent = []
            Ns = []
            println("Computing entanglement for nodes at distance $(distance-1), λ = $λmean")
            println()
            println()
            flush(stdout)

            for (v1, v2) in pairs[distance]
                sol = find_path(g, v1, v2, samples=256, σ=0.15, λmean=λmean, progressbar=false)
                if !isnothing(sol)
                    path, λ, N = sol
                    println("Path to connect $v1 to $v2:")
                    println(path)
                    println("Entanglement: $(2*(1-λ))")
                    println()
                    println()
                    flush(stdout)
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
            write(file, "$(dist-1) $(mean_ent[dist]) $(mean_N[dist]) \n")
        end
        flush(file)
        close(file)
    end
end

function run_dev(n, n_samples, σ)

    if !isdir("out_files/$(n)")
        mkdir("out_files/$(n)")
    end
    s_num = isempty(readdir("out_files/$n")) ? 1 : (maximum([parse(Int64, sample) for sample in readdir("out_files/$n")]) + 1)
    for i in s_num:s_num+n_samples-1
        if !isdir("out_files/$n/$i")
            mkdir("out_files/$n/$i")
        end
    end

    λs = collect(range(0.5, stop=1.0, length=101))
    for j in 1:n_samples
        for (i, λmean) in enumerate(λs)
            g = noisy_square_lattice(n, σ=σ, λmean=λmean)

            pairs = pairs_by_distance(g)

            λmean = mean(get_weights(g, pair...) for pair in pairs[1])

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

# run_mean(6)
# run_dev(6, 1, 0.01)