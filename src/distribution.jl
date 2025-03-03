using Random
using DataStructures
using GraphPlot
using GraphRecipes
using Plots
using StatsPlots
using Trapz
using Distributions
using Statistics
using Combinatorics
using LinearAlgebra
using FileIO
using ImageIO
using Roots

# P(λ, δ) = 2*(2*(1 - λ))^(1/δ - 1)/δ

# sample_λ(δ) = round(0.5*(rand(Beta(1, 1/δ)) + 1), digits=2)

function sample_λ(N, x_min, x_max, x_mean, std_dev)
    if !(x_min <= x_mean <= x_max)
        error("x_mean must be within the range [x_min, x_max]")
    end

    if x_mean == x_min || x_mean == x_max
        return fill(x_mean, N)
    end

    # Function to solve for the true underlying normal mean μ_true
    function find_true_mean(mu)
        truncated_dist = Truncated(Normal(mu, std_dev), x_min, x_max)
        return mean(truncated_dist) - x_mean
    end

    # Use a numerical solver to find μ_true such that the truncated mean equals x_mean
    μ_true = find_zero(find_true_mean, x_mean)

    # Define the correctly parameterized truncated normal distribution
    truncated_dist = Truncated(Normal(μ_true, std_dev), x_min, x_max)

    # Sample N values from this distribution
    return rand(truncated_dist, N)
end

function test_distribution()
    N = 10000
    x_min, x_max = 0.5, 1
    std_dev_1, std_dev_2 = 0.07, 0.01
    mean_range = range(0.5, stop=1, length=101)[2:end-1]  # Interpolate means between 0.5 and 1

    # Initialize the plot with FFMPEG backend for video creation
    output_video = "distribution_evolution_ffmpeg.mp4"
    fps = 5  # Frames per second
    anim = @animate for x_mean in mean_range
        values1 = sample_λ(N, x_min, x_max, x_mean, std_dev_1)
        println("Predicted mean: ", x_mean)
        println("Mean of values 1: ", mean(values1))

        values2 = sample_λ(N, x_min, x_max, x_mean, std_dev_2)
        println("Mean of values 2: ", mean(values2))
        
        # Plot the PDFs
        p = plot(legend=:topright, xlabel="λ", ylabel="Density", title="PDF of Gaussian Distributions for Quenches")
        density!(values1, label="σ = 0.07", color=:skyblue, lw=2, fill=(0,:skyblue,0.4))
        density!(values2, label="σ = 0.01", color=:lightgreen, lw=2, fill=(0,:lightgreen,0.4))
        vline!([x_mean], linestyle=:dash, color=:red, label="Mean: $x_mean")

        xlims!(x_min, x_max)
        ylims!(0, 50)

        println()
        println()
    end

    # Save the animation as a video using FFMPEG
    mp4(anim, output_video, fps=fps)

    println("Video created: $output_video")
end

test_distribution()