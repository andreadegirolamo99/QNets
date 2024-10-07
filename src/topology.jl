include("graph.jl")
include("distribution.jl")


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

function noisy_square_lattice(n::Int64; λmean=0.5, σ=0.01)
    g = Graph(n^2)

    λs = generate_gaussian_values(2*n*(n-1), 0.5, 1, λmean, σ)

    for j in 1:n
        for i in 1:n-1
            add_edge!(g, i+n*(j-1), i+n*(j-1)+1, splice!(λs, 1:length(λs)))
        end
    end
    for j in 1:n-1
        for i in 1:n
            add_edge!(g, i+n*(j-1), i+n*j, splice!(λs, 1:length(λs)))
        end
    end

    return g
end

function noisy_diagonal_square_lattice(n::Int64; λmean=0.5, σ=0.01)
    g = Graph(n^2)

    λs = generate_gaussian_values(4*n^2 - 6*n + 2, 0.5, 1, λmean, σ)

    for j in 1:n-1
        for i in 1:n-1
            add_edge!(g, i+n*(j-1), i+n*(j-1)+1, splice!(λs, 1:length(λs)))
            add_edge!(g, i+n*(j-1), i+n*j+1, splice!(λs, 1:length(λs)))
        end
    end
    for i in 1:n-1
        add_edge!(g, i+n*(n-1), i+n*(n-1)+1, splice!(λs, 1:length(λs)))
    end

    for j in 1:n-1
        for i in 1:n-1
            add_edge!(g, i+n*(j-1), i+n*j, splice!(λs, 1:length(λs)))
            add_edge!(g, i+n*(j-1)+1, i+n*j, splice!(λs, 1:length(λs)))
        end
        add_edge!(g, n*j, n*(j+1), splice!(λs, 1:length(λs)))
    end

    return g
end