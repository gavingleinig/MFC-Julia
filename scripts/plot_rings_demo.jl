using Random
using Plots
using MetricForestCompletion

include("plot_utils_convex_hull.jl")  # for to_2d, used only if dim > 2

function plot_rings(points::Vector{Vector{T}}, labels::Vector{Int}; title::String="Generated Rings") where {T}
    dim = length(points[1])
    coords = dim == 2 ? points : to_2d(points)

    xs = [c[1] for c in coords]
    ys = [c[2] for c in coords]
    group_labels = string.("Ring ", labels)

    scatter(xs, ys; group = group_labels,
            xlabel = dim == 2 ? "x1" : "PC1",
            ylabel = dim == 2 ? "x2" : "PC2",
            title = title, markersize = 3, markerstrokewidth = 0,
            aspect_ratio = :equal, legend = :outertopright)
end

function main()
    rng = Xoshiro(123)

    dim             = 2
    num_rings       = 6
    points_per_ring = 333
    radius_range    = (2.0, 10.0)
    noise_sigma     = 0.3   # try 0.0 for perfectly thin rings, or bump this up
                             # to see rings start to overlap/blur together

    points, labels = generate_rings(
        Float64;
        dim = dim,
        num_rings = num_rings,
        points_per_ring = points_per_ring,
        radius_range = radius_range,
        noise_sigma = noise_sigma,
        rng = rng,
    )

    plt = plot_rings(points, labels; title = "Concentric rings (noise_sigma=$noise_sigma)")

    outpath = joinpath(@__DIR__, "rings_demo.png")
    savefig(plt, outpath)
    println("Saved: ", outpath)
    return plt
end

main()
