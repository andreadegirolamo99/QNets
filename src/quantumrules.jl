include("graph.jl")


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
    α, β = sort(get_weights(g, v1, v2))[1:2]

    λ = λdistill(α, β)

    rem_edge!(g, min(v1, v2), max(v1, v2), α)
    rem_edge!(g, min(v1, v2), max(v1, v2), β)
    add_edge!(g, min(v1, v2), max(v1, v2), λ)
end

function connect_ends(g::QGraph)
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

function apply(g::QGraph, path)
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