using Random
using Distributions
using LinearAlgebra
using Plots
using MetricForestCompletion

#include("generate_gaussians.jl")        # fixed version -> generate_gaussians(...)
#include("generate_gaussians_buggy.jl")  # original version -> generate_gaussians_buggy(...)

function scatter_2d(points::Vector{Vector{T}}, labels::Vector{Int}; title::String) where {T}
    xs = [p[1] for p in points]
    ys = [p[2] for p in points]
    group_labels = string.("Cluster ", labels)
    scatter(xs, ys; group=group_labels,
            xlabel="x1", ylabel="x2", title=title,
            markersize=3, markerstrokewidth=0, legend=:outertopright,
            aspect_ratio=:equal)
end

function main()
    dim               = 2          # kept at 2 so the diagonal collapse is directly visible
    num_centers       = 5
    points_per_center = 200
    mean_range        = (-10.0, 10.0)
    sigma_range       = (0.5, 2.0)

    # --- Original / buggy generator ---
    buggy_points, buggy_labels = generate_gaussians(
        Float64;
        dim = dim,
        num_centers = num_centers,
        points_per_center = points_per_center,
        mean_range = mean_range,
        sigma_range = sigma_range,
    )

    # --- Fixed generator (for side-by-side comparison) ---
    # fixed_points, fixed_labels = generate_gaussians(
    #     Float64;
    #     dim = dim,
    #     num_centers = num_centers,
    #     points_per_center = points_per_center,
    #     mean_range = mean_range,
    #     sigma_range = sigma_range,
    #     rng = Xoshiro(123),
    # )

    plt_buggy = scatter_2d(buggy_points, buggy_labels; title = "Original: centers collapse near x1 = x2")
    # plt_fixed = scatter_2d(fixed_points, fixed_labels; title = "Fixed: centers scattered in the plane")

    combined = plot(plt_buggy)

    outpath = joinpath(@__DIR__, "gaussian_bug_comparison.png")
    savefig(combined, outpath)
    println("Saved comparison plot to: ", outpath)

    return combined
end

main()