function square_lattice(n::Int64; λ=0.5)
    g = QGraph(n^2)

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
    g = QGraph(n^2)

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

function honeycomb(n::Int64; λ=0.5)
    g = QGraph(2*n*(n+2))

    prev_hex_vs = []
    for i in 1:(div(n, 2) + mod(n, 2))
        add_edge!(g, 6*i - 5, 6*i - 4, λ)
        add_edge!(g, 6*i - 5, 6*i - 3, λ)
        add_edge!(g, 6*i - 4, 6*i - 2, λ)
        add_edge!(g, 6*i - 3, 6*i - 1, λ)
        add_edge!(g, 6*i - 2, 6*i, λ)
        add_edge!(g, 6*i - 1, 6*i, λ)
        push!(prev_hex_vs, [6*i-1, 6*i])
    end

    for i in 1:(div(n, 2) + mod(n, 2) - 1)
        add_edge!(g, 6*i - 2, 6*i + 3, λ)
    end

    ends = []

    if mod(n, 2) == 0
        add_edge!(g, 6*(div(n, 2) + mod(n, 2)) - 2, 6*(div(n, 2) + mod(n, 2)) + 1, λ)
        push!(ends, 6*(div(n, 2) + mod(n, 2)) + 1)
    end

    vertices = collect((6*(div(n, 2) + mod(n, 2)) + 2 - mod(n, 2)):2*n*(n+2))

    final_vs = []
    vertex_end = nothing

    for j in 1:n-1
        v_useful = nothing
        v_intermediate = []
        for i in 1:(div(n, 2) + mod(n, 2))
            v1, v2, v3, v4 = splice!(vertices, 1:4)
            if n > 2
                if i == 1
                    push!(v_intermediate, [])
                    push!(v_intermediate[end], v2)
                elseif i == (div(n, 2) + mod(n, 2))
                    push!(v_intermediate[end], v1)
                else
                    push!(v_intermediate[end], v1)
                    push!(v_intermediate, [])
                    push!(v_intermediate[end], v2)
                end
            end
            if i == (div(n, 2) + mod(n, 2))
                v_useful = v2
            end
            add_edge!(g, prev_hex_vs[i][1], v1, λ)
            add_edge!(g, v1, v3, λ)
            add_edge!(g, v3, v4, λ)
            add_edge!(g, prev_hex_vs[i][2], v2, λ)
            add_edge!(g, v2, v4, λ)
            prev_hex_vs[i] = [v3, v4]

            if j == n-1
                if (div(n, 2) + mod(n, 2)) == 1
                    vertex_end = v4
                elseif i == 1
                    push!(final_vs, [])
                    push!(final_vs[end], v4)
                elseif i == (div(n, 2) + mod(n, 2))
                    push!(final_vs[end], v3)
                    vertex_end = v4
                else
                    push!(final_vs[end], v3)
                    push!(final_vs, [])
                    push!(final_vs[end], v4)
                end
            end
        end

        for i in eachindex(v_intermediate)
            add_edge!(g, v_intermediate[i][1], v_intermediate[i][2], λ)
        end

        if mod(n, 2) == 0
            v = splice!(vertices, 1)
            add_edge!(g, v_useful, v, λ)
            push!(ends, v)
        end
    end

    for i in 1:(length(ends)-1)
        v = splice!(vertices, 1)
        add_edge!(g, ends[i], v, λ)
        add_edge!(g, ends[i+1], v, λ)
    end

    for i in eachindex(final_vs)
        v1, v2 = splice!(vertices, 1:2)
        add_edge!(g, final_vs[i][1], v1, λ)
        add_edge!(g, v1, v2, λ)
        add_edge!(g, v2, final_vs[i][2], λ)
    end

    if mod(n, 2) == 0
        v1, v2, v3 = splice!(vertices, 1:3)
        add_edge!(g, vertex_end, v1, λ)
        add_edge!(g, v1, v2, λ)
        add_edge!(g, v2, v3, λ)
        add_edge!(g, v3, ends[end], λ)
    end

    return g
end

function fully_connected_honeycomb(n::Int64; λ=0.5)
    g = QGraph(2*n*(n+2))

    prev_hex_vs = []

    hexagon = [ [0, 0, 0, 0, 0, 0] for i in 1:n, j in 1:n ]
    for i in 1:(div(n, 2) + mod(n, 2))
        hexagon[1,2*(i-1) + 1][1] = 6*i - 5
        hexagon[1,2*(i-1) + 1][2] = 6*i - 4
        hexagon[1,2*(i-1) + 1][3] = 6*i - 3
        hexagon[1,2*(i-1) + 1][4] = 6*i - 2
        hexagon[1,2*(i-1) + 1][5] = 6*i - 1
        hexagon[1,2*(i-1) + 1][6] = 6*i
        if size(hexagon)[1] > 1
            hexagon[2, 2*(i-1) + 1][1] = 6*i - 1
            hexagon[2, 2*(i-1) + 1][2] = 6*i
        end
        if size(hexagon)[2] > 2*(i-1) + 1
            hexagon[1,2*(i-1) + 2][1] = 6*i - 2
            hexagon[1,2*(i-1) + 2][3] = 6*i
        end
        if i > 1
            hexagon[1,2*(i-1)][2] = 6*i - 3
            hexagon[1,2*(i-1)][4] = 6*i - 1
        end

        add_edge!(g, 6*i - 5, 6*i - 4, λ)
        add_edge!(g, 6*i - 5, 6*i - 3, λ)
        add_edge!(g, 6*i - 4, 6*i - 2, λ)
        add_edge!(g, 6*i - 3, 6*i - 1, λ)
        add_edge!(g, 6*i - 2, 6*i, λ)
        add_edge!(g, 6*i - 1, 6*i, λ)

        push!(prev_hex_vs, [6*i-1, 6*i])
    end

    for i in 1:(div(n, 2) + mod(n, 2) - 1)
        add_edge!(g, 6*i - 2, 6*i + 3, λ)
    end

    ends = []

    if mod(n, 2) == 0
        add_edge!(g, 6*(div(n, 2) + mod(n, 2)) - 2, 6*(div(n, 2) + mod(n, 2)) + 1, λ)
        hexagon[1, size(hexagon)[2]][2] = 6*(div(n, 2) + mod(n, 2)) + 1
        push!(ends, 6*(div(n, 2) + mod(n, 2)) + 1)
    end

    vertices = collect((6*(div(n, 2) + mod(n, 2)) + 2 - mod(n, 2)):2*n*(n+2))

    final_vs = []
    vertex_end = nothing

    for j in 1:n-1
        v_useful = nothing
        v_intermediate = []
        for i in 1:(div(n, 2) + mod(n, 2))
            v1, v2, v3, v4 = splice!(vertices, 1:4)
            hexagon[j+1, 2*(i-1)+1][3] = v1
            hexagon[j+1, 2*(i-1)+1][4] = v2
            hexagon[j+1, 2*(i-1)+1][5] = v3
            hexagon[j+1, 2*(i-1)+1][6] = v4


            if size(hexagon)[1] > j+1
                hexagon[j+2, 2*(i-1)+1][1] = v3
                hexagon[j+2, 2*(i-1)+1][2] = v4
            end
            if size(hexagon)[2] > 2*(i-1) + 1
                hexagon[j+1, 2*(i-1) + 2][1] = v2
                hexagon[j+1, 2*(i-1) + 2][3] = v4
                hexagon[j, 2*(i-1) + 2][5] = v2
            end
            if i > 1
                hexagon[j+1,2*(i-1)][2] = v1
                hexagon[j,2*(i-1)][6] = v1
                hexagon[j+1,2*(i-1)][4] = v3
            end

            if n > 2
                if i == 1
                    push!(v_intermediate, [])
                    push!(v_intermediate[end], v2)
                elseif i == (div(n, 2) + mod(n, 2))
                    push!(v_intermediate[end], v1)
                else
                    push!(v_intermediate[end], v1)
                    push!(v_intermediate, [])
                    push!(v_intermediate[end], v2)
                end
            end
            if i == (div(n, 2) + mod(n, 2))
                v_useful = v2
            end
            add_edge!(g, prev_hex_vs[i][1], v1, λ)
            add_edge!(g, v1, v3, λ)
            add_edge!(g, v3, v4, λ)
            add_edge!(g, prev_hex_vs[i][2], v2, λ)
            add_edge!(g, v2, v4, λ)
            prev_hex_vs[i] = [v3, v4]

            if j == n-1
                if (div(n, 2) + mod(n, 2)) == 1
                    vertex_end = v4
                elseif i == 1
                    push!(final_vs, [])
                    push!(final_vs[end], v4)
                elseif i == (div(n, 2) + mod(n, 2))
                    push!(final_vs[end], v3)
                    vertex_end = v4
                else
                    push!(final_vs[end], v3)
                    push!(final_vs, [])
                    push!(final_vs[end], v4)
                end
            end
        end

        for i in eachindex(v_intermediate)
            add_edge!(g, v_intermediate[i][1], v_intermediate[i][2], λ)
        end

        if mod(n, 2) == 0
            v = splice!(vertices, 1)
            add_edge!(g, v_useful, v, λ)
            push!(ends, v)
            hexagon[j+1, size(hexagon)[2]][2] = v
            hexagon[j, size(hexagon)[2]][6] = v
        end
    end

    for i in 1:(length(ends)-1)
        v = splice!(vertices, 1)
        add_edge!(g, ends[i], v, λ)
        add_edge!(g, ends[i+1], v, λ)
        hexagon[i, size(hexagon)[2]][4] = v
    end

    for i in eachindex(final_vs)
        v1, v2 = splice!(vertices, 1:2)
        hexagon[size(hexagon)[1], 2*i][5] = v1
        hexagon[size(hexagon)[1], 2*i][6] = v2
        add_edge!(g, final_vs[i][1], v1, λ)
        add_edge!(g, v1, v2, λ)
        add_edge!(g, v2, final_vs[i][2], λ)
    end

    if mod(n, 2) == 0
        v1, v2, v3 = splice!(vertices, 1:3)
        add_edge!(g, v2, v3, λ)
        hexagon[size(hexagon)[1], size(hexagon)[2]][4] = v1
        hexagon[size(hexagon)[1], size(hexagon)[2]][5] = v2
        hexagon[size(hexagon)[1], size(hexagon)[2]][6] = v3
    end

    for hex in hexagon
        add_edge!(g, hex[1], hex[4], λ)
        add_edge!(g, hex[1], hex[5], λ)
        add_edge!(g, hex[1], hex[6], λ)
        add_edge!(g, hex[2], hex[3], λ)
        add_edge!(g, hex[2], hex[5], λ)
        add_edge!(g, hex[2], hex[6], λ)
        add_edge!(g, hex[3], hex[4], λ)
        add_edge!(g, hex[3], hex[6], λ)
        add_edge!(g, hex[4], hex[5], λ)
    end

    return g
end

function noisy_square_lattice(n::Int64; λmean=0.5, σ=0.01)
    g = QGraph(n^2)

    λs = sample_λ(2*n*(n-1), 0.5, 1, λmean, σ)

    for j in 1:n
        for i in 1:n-1
            add_edge!(g, i+n*(j-1), i+n*(j-1)+1, splice!(λs, rand(1:length(λs))))
        end
    end
    for j in 1:n-1
        for i in 1:n
            add_edge!(g, i+n*(j-1), i+n*j, splice!(λs, rand(1:length(λs))))
        end
    end

    return g
end

function noisy_diagonal_square_lattice(n::Int64; λmean=0.5, σ=0.01)
    g = QGraph(n^2)

    λs = sample_λ(4*n^2 - 6*n + 2, 0.5, 1, λmean, σ)

    for j in 1:n-1
        for i in 1:n-1
            add_edge!(g, i+n*(j-1), i+n*(j-1)+1, splice!(λs, rand(1:length(λs))))
            add_edge!(g, i+n*(j-1), i+n*j+1, splice!(λs, rand(1:length(λs))))
        end
    end
    for i in 1:n-1
        add_edge!(g, i+n*(n-1), i+n*(n-1)+1, splice!(λs, rand(1:length(λs))))
    end

    for j in 1:n-1
        for i in 1:n-1
            add_edge!(g, i+n*(j-1), i+n*j, splice!(λs, rand(1:length(λs))))
            add_edge!(g, i+n*(j-1)+1, i+n*j, splice!(λs, rand(1:length(λs))))
        end
        add_edge!(g, n*j, n*(j+1), splice!(λs, rand(1:length(λs))))
    end

    return g
end

function noisy_honeycomb(n::Int64; λmean=0.5, σ=0.01)
    g = QGraph(2*n*(n+2))

    λs = sample_λ(3*n^2 + 4*n - 1, 0.5, 1, λmean, σ)

    prev_hex_vs = []
    for i in 1:(div(n, 2) + mod(n, 2))
        add_edge!(g, 6*i - 5, 6*i - 4, splice!(λs, rand(1:length(λs))))
        add_edge!(g, 6*i - 5, 6*i - 3, splice!(λs, rand(1:length(λs))))
        add_edge!(g, 6*i - 4, 6*i - 2, splice!(λs, rand(1:length(λs))))
        add_edge!(g, 6*i - 3, 6*i - 1, splice!(λs, rand(1:length(λs))))
        add_edge!(g, 6*i - 2, 6*i, splice!(λs, rand(1:length(λs))))
        add_edge!(g, 6*i - 1, 6*i, splice!(λs, rand(1:length(λs))))
        push!(prev_hex_vs, [6*i-1, 6*i])
    end

    for i in 1:(div(n, 2) + mod(n, 2) - 1)
        add_edge!(g, 6*i - 2, 6*i + 3, splice!(λs, rand(1:length(λs))))
    end

    ends = []

    if mod(n, 2) == 0
        add_edge!(g, 6*(div(n, 2) + mod(n, 2)) - 2, 6*(div(n, 2) + mod(n, 2)) + 1, splice!(λs, rand(1:length(λs))))
        push!(ends, 6*(div(n, 2) + mod(n, 2)) + 1)
    end

    vertices = collect((6*(div(n, 2) + mod(n, 2)) + 2 - mod(n, 2)):2*n*(n+2))

    final_vs = []
    vertex_end = nothing

    for j in 1:n-1
        v_useful = nothing
        v_intermediate = []
        for i in 1:(div(n, 2) + mod(n, 2))
            v1, v2, v3, v4 = splice!(vertices, 1:4)
            if n > 2
                if i == 1
                    push!(v_intermediate, [])
                    push!(v_intermediate[end], v2)
                elseif i == (div(n, 2) + mod(n, 2))
                    push!(v_intermediate[end], v1)
                else
                    push!(v_intermediate[end], v1)
                    push!(v_intermediate, [])
                    push!(v_intermediate[end], v2)
                end
            end
            if i == (div(n, 2) + mod(n, 2))
                v_useful = v2
            end
            add_edge!(g, prev_hex_vs[i][1], v1, splice!(λs, rand(1:length(λs))))
            add_edge!(g, v1, v3, splice!(λs, rand(1:length(λs))))
            add_edge!(g, v3, v4, splice!(λs, rand(1:length(λs))))
            add_edge!(g, prev_hex_vs[i][2], v2, splice!(λs, rand(1:length(λs))))
            add_edge!(g, v2, v4, splice!(λs, rand(1:length(λs))))
            prev_hex_vs[i] = [v3, v4]

            if j == n-1
                if (div(n, 2) + mod(n, 2)) == 1
                    vertex_end = v4
                elseif i == 1
                    push!(final_vs, [])
                    push!(final_vs[end], v4)
                elseif i == (div(n, 2) + mod(n, 2))
                    push!(final_vs[end], v3)
                    vertex_end = v4
                else
                    push!(final_vs[end], v3)
                    push!(final_vs, [])
                    push!(final_vs[end], v4)
                end
            end
        end

        for i in eachindex(v_intermediate)
            add_edge!(g, v_intermediate[i][1], v_intermediate[i][2], splice!(λs, rand(1:length(λs))))
        end

        if mod(n, 2) == 0
            v = splice!(vertices, 1)
            add_edge!(g, v_useful, v, splice!(λs, rand(1:length(λs))))
            push!(ends, v)
        end
    end

    for i in 1:(length(ends)-1)
        v = splice!(vertices, 1)
        add_edge!(g, ends[i], v, splice!(λs, rand(1:length(λs))))
        add_edge!(g, ends[i+1], v, splice!(λs, rand(1:length(λs))))
    end

    for i in eachindex(final_vs)
        v1, v2 = splice!(vertices, 1:2)
        add_edge!(g, final_vs[i][1], v1, splice!(λs, rand(1:length(λs))))
        add_edge!(g, v1, v2, splice!(λs, rand(1:length(λs))))
        add_edge!(g, v2, final_vs[i][2], splice!(λs, rand(1:length(λs))))
    end

    if mod(n, 2) == 0
        v1, v2, v3 = splice!(vertices, 1:3)
        add_edge!(g, vertex_end, v1, splice!(λs, rand(1:length(λs))))
        add_edge!(g, v1, v2, splice!(λs, rand(1:length(λs))))
        add_edge!(g, v2, v3, splice!(λs, rand(1:length(λs))))
        add_edge!(g, v3, ends[end], splice!(λs, rand(1:length(λs))))
    end

    return g
end

function noisy_fully_connected_honeycomb(n::Int64; λmean=0.5, σ=0.01)
    g = QGraph(2*n*(n+2))

    λs = sample_λ(12*n^2 + 4*n - 1, 0.5, 1, λmean, σ)

    prev_hex_vs = []

    hexagon = [ [0, 0, 0, 0, 0, 0] for i in 1:n, j in 1:n ]
    for i in 1:(div(n, 2) + mod(n, 2))
        hexagon[1,2*(i-1) + 1][1] = 6*i - 5
        hexagon[1,2*(i-1) + 1][2] = 6*i - 4
        hexagon[1,2*(i-1) + 1][3] = 6*i - 3
        hexagon[1,2*(i-1) + 1][4] = 6*i - 2
        hexagon[1,2*(i-1) + 1][5] = 6*i - 1
        hexagon[1,2*(i-1) + 1][6] = 6*i
        if size(hexagon)[1] > 1
            hexagon[2, 2*(i-1) + 1][1] = 6*i - 1
            hexagon[2, 2*(i-1) + 1][2] = 6*i
        end
        if size(hexagon)[2] > 2*(i-1) + 1
            hexagon[1,2*(i-1) + 2][1] = 6*i - 2
            hexagon[1,2*(i-1) + 2][3] = 6*i
        end
        if i > 1
            hexagon[1,2*(i-1)][2] = 6*i - 3
            hexagon[1,2*(i-1)][4] = 6*i - 1
        end

        add_edge!(g, 6*i - 5, 6*i - 4, splice!(λs, rand(1:length(λs))))
        add_edge!(g, 6*i - 5, 6*i - 3, splice!(λs, rand(1:length(λs))))
        add_edge!(g, 6*i - 4, 6*i - 2, splice!(λs, rand(1:length(λs))))
        add_edge!(g, 6*i - 3, 6*i - 1, splice!(λs, rand(1:length(λs))))
        add_edge!(g, 6*i - 2, 6*i, splice!(λs, rand(1:length(λs))))
        add_edge!(g, 6*i - 1, 6*i, splice!(λs, rand(1:length(λs))))

        push!(prev_hex_vs, [6*i-1, 6*i])
    end

    for i in 1:(div(n, 2) + mod(n, 2) - 1)
        add_edge!(g, 6*i - 2, 6*i + 3, splice!(λs, rand(1:length(λs))))
    end

    ends = []

    if mod(n, 2) == 0
        add_edge!(g, 6*(div(n, 2) + mod(n, 2)) - 2, 6*(div(n, 2) + mod(n, 2)) + 1, splice!(λs, rand(1:length(λs))))
        hexagon[1, size(hexagon)[2]][2] = 6*(div(n, 2) + mod(n, 2)) + 1
        push!(ends, 6*(div(n, 2) + mod(n, 2)) + 1)
    end

    vertices = collect((6*(div(n, 2) + mod(n, 2)) + 2 - mod(n, 2)):2*n*(n+2))

    final_vs = []
    vertex_end = nothing

    for j in 1:n-1
        v_useful = nothing
        v_intermediate = []
        for i in 1:(div(n, 2) + mod(n, 2))
            v1, v2, v3, v4 = splice!(vertices, 1:4)
            hexagon[j+1, 2*(i-1)+1][3] = v1
            hexagon[j+1, 2*(i-1)+1][4] = v2
            hexagon[j+1, 2*(i-1)+1][5] = v3
            hexagon[j+1, 2*(i-1)+1][6] = v4


            if size(hexagon)[1] > j+1
                hexagon[j+2, 2*(i-1)+1][1] = v3
                hexagon[j+2, 2*(i-1)+1][2] = v4
            end
            if size(hexagon)[2] > 2*(i-1) + 1
                hexagon[j+1, 2*(i-1) + 2][1] = v2
                hexagon[j+1, 2*(i-1) + 2][3] = v4
                hexagon[j, 2*(i-1) + 2][5] = v2
            end
            if i > 1
                hexagon[j+1,2*(i-1)][2] = v1
                hexagon[j,2*(i-1)][6] = v1
                hexagon[j+1,2*(i-1)][4] = v3
            end

            if n > 2
                if i == 1
                    push!(v_intermediate, [])
                    push!(v_intermediate[end], v2)
                elseif i == (div(n, 2) + mod(n, 2))
                    push!(v_intermediate[end], v1)
                else
                    push!(v_intermediate[end], v1)
                    push!(v_intermediate, [])
                    push!(v_intermediate[end], v2)
                end
            end
            if i == (div(n, 2) + mod(n, 2))
                v_useful = v2
            end
            add_edge!(g, prev_hex_vs[i][1], v1, splice!(λs, rand(1:length(λs))))
            add_edge!(g, v1, v3, splice!(λs, rand(1:length(λs))))
            add_edge!(g, v3, v4, splice!(λs, rand(1:length(λs))))
            add_edge!(g, prev_hex_vs[i][2], v2, splice!(λs, rand(1:length(λs))))
            add_edge!(g, v2, v4, splice!(λs, rand(1:length(λs))))
            prev_hex_vs[i] = [v3, v4]

            if j == n-1
                if (div(n, 2) + mod(n, 2)) == 1
                    vertex_end = v4
                elseif i == 1
                    push!(final_vs, [])
                    push!(final_vs[end], v4)
                elseif i == (div(n, 2) + mod(n, 2))
                    push!(final_vs[end], v3)
                    vertex_end = v4
                else
                    push!(final_vs[end], v3)
                    push!(final_vs, [])
                    push!(final_vs[end], v4)
                end
            end
        end

        for i in eachindex(v_intermediate)
            add_edge!(g, v_intermediate[i][1], v_intermediate[i][2], splice!(λs, rand(1:length(λs))))
        end

        if mod(n, 2) == 0
            v = splice!(vertices, 1)
            add_edge!(g, v_useful, v, splice!(λs, rand(1:length(λs))))
            push!(ends, v)
            hexagon[j+1, size(hexagon)[2]][2] = v
            hexagon[j, size(hexagon)[2]][6] = v
        end
    end

    for i in 1:(length(ends)-1)
        v = splice!(vertices, 1)
        add_edge!(g, ends[i], v, splice!(λs, rand(1:length(λs))))
        add_edge!(g, ends[i+1], v, splice!(λs, rand(1:length(λs))))
        hexagon[i, size(hexagon)[2]][4] = v
    end

    for i in eachindex(final_vs)
        v1, v2 = splice!(vertices, 1:2)
        hexagon[size(hexagon)[1], 2*i][5] = v1
        hexagon[size(hexagon)[1], 2*i][6] = v2
        add_edge!(g, final_vs[i][1], v1, splice!(λs, rand(1:length(λs))))
        add_edge!(g, v1, v2, splice!(λs, rand(1:length(λs))))
        add_edge!(g, v2, final_vs[i][2], splice!(λs, rand(1:length(λs))))
    end

    if mod(n, 2) == 0
        v1, v2, v3 = splice!(vertices, 1:3)
        add_edge!(g, v2, v3, splice!(λs, rand(1:length(λs))))
        hexagon[size(hexagon)[1], size(hexagon)[2]][4] = v1
        hexagon[size(hexagon)[1], size(hexagon)[2]][5] = v2
        hexagon[size(hexagon)[1], size(hexagon)[2]][6] = v3
    end

    for hex in hexagon
        add_edge!(g, hex[1], hex[4], splice!(λs, rand(1:length(λs))))
        add_edge!(g, hex[1], hex[5], splice!(λs, rand(1:length(λs))))
        add_edge!(g, hex[1], hex[6], splice!(λs, rand(1:length(λs))))
        add_edge!(g, hex[2], hex[3], splice!(λs, rand(1:length(λs))))
        add_edge!(g, hex[2], hex[5], splice!(λs, rand(1:length(λs))))
        add_edge!(g, hex[2], hex[6], splice!(λs, rand(1:length(λs))))
        add_edge!(g, hex[3], hex[4], splice!(λs, rand(1:length(λs))))
        add_edge!(g, hex[3], hex[6], splice!(λs, rand(1:length(λs))))
        add_edge!(g, hex[4], hex[5], splice!(λs, rand(1:length(λs))))
    end

    return g
end