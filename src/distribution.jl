using Random
using Distributions
using Statistics
using Combinatorics
using LinearAlgebra
using FileIO
using ImageIO

# P(λ, δ) = 2*(2*(1 - λ))^(1/δ - 1)/δ

# sample_λ(δ) = round(0.5*(rand(Beta(1, 1/δ)) + 1), digits=2)

function sample_λ(N, x_min, x_max, x_mean, std_dev)
    if !(x_min <= x_mean <= x_max)
        error("x_mean must be within the range [x_min, x_max]")
    end

    if x_mean == x_min || x_mean == x_max
        return fill(x_mean, N)
    end
    
    # Step 1: Generate N random values from a Gaussian distribution
    values = rand(Normal(x_mean, std_dev), N)
    
    # Step 2: Clip values to ensure they remain within the range [x_min, x_max]
    clipped_values = clamp.(values, x_min, x_max)
    
    # Step 3: Calculate the current mean of the clipped values
    current_mean = mean(clipped_values)
    
    # Step 4: Adjust the clipped values to ensure the mean is exactly x_mean
    diff = x_mean - current_mean
    adjusted_values = clipped_values .+ diff
    
    # Ensure after adjustment that values don't exceed bounds
    adjusted_values = clamp.(adjusted_values, x_min, x_max)
    
    # Step 5: Final check and adjustment if the mean is still not x_mean
    final_mean = mean(adjusted_values)
    max_iters = 10000
    cnt = 0
    while abs(final_mean - x_mean) > 1e-16 && cnt <= max_iters
        # Scale the difference proportionally to the range
        diff = x_mean - final_mean
        adjusted_values .= adjusted_values .+ diff
        
        # Clip again to avoid out of range values
        adjusted_values = clamp.(adjusted_values, x_min, x_max)
        
        # Recompute the mean
        final_mean = mean(adjusted_values)
        cnt += 1
    end
    
    return adjusted_values
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
        values1 = generate_gaussian_values(N, x_min, x_max, x_mean, std_dev_1)
        println("Predicted mean: ", x_mean)
        println("Mean of values 1: ", mean(values1))

        values2 = generate_gaussian_values(N, x_min, x_max, x_mean, std_dev_2)
        println("Mean of values 2: ", mean(values2))
        
        # Plot the PDFs
        p = plot(legend=:topright, xlabel="λ", ylabel="Density", title="PDF of Gaussian Distributions for Quenches")
        density!(values1, label="σ = 0.07", color=:skyblue, lw=2, fill=(0,:skyblue,0.4))
        density!(values2, label="σ = 0.01", color=:lightgreen, lw=2, fill=(0,:lightgreen,0.4))
        vline!([x_mean], linestyle=:dash, color=:red, label="Mean: $x_mean")
        
        xlims!(x_min, x_max)
        ylims!(0, 50)

        kde1 = density(values1)
        kde2 = density(values2)

        integral1 = trapz(kde1.x, kde1.y)
        integral2 = trapz(kde2.x, kde2.y)
        println()

        println("Integral of distribution 1: $integral1")
        println("Integral of distribution 2: $integral2")

        println()
        println()
    end

    # Save the animation as a video using FFMPEG
    mp4(anim, output_video, fps=fps)

    println("Video created: $output_video")
end