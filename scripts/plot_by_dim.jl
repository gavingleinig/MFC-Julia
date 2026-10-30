using CSV
using DataFrames
using Plots
using Statistics

# Accuracy vs. dimension: one panel per Gaussian (or ring) count, one line per
# method, shared legend underneath. Reads the same CSVs as plot_grid.jl.
#
# `sigma` picks which noise level to show: the row whose Sigma is closest to it.
# With sigma = nothing, each point is the mean over every sigma in the sweep.
#  gaussians = [3, 4, 5]
#     # gaussians = [16, 32, 64, 128, 256]  # for full sweep
 
#     # change number of dimensions
#     # dims = [2, 3, 4]
#     dims = [4, 16, 64, 128]  # for full sweep


function plot_by_dimension(;
    dims::Vector{Int} = [4, 16, 64],
    gaussians::Vector{Int} = [16, 64],
    metric::String = "NMI",                     # "ARI" or "NMI"
    sigma::Union{Real, Nothing} = nothing,
    input_dir::String = "data/Systhetic/",
    name_prefix::String = "final_clustering_results",
    name_infix::String = "gaus",                # "gaus" or "rings" in the file name
    count_label::String = "Gaussians",          # panel title label
    output_path::String = joinpath(input_dir, "by_dim_$(metric).png"),
)
    methods = [
        ("KC",          "K-Centering", :circle),
        ("MST",         "Optimal MST", :square),
        ("Naive",       "Naive ST",    :dtriangle),
        ("MFC_Approx",  "MFC Approx",  :diamond),
        ("MFC_Optimal", "MFC Optimal", :hexagon),
        ("MFC_Simple",  "MFC Simple",  :xcross),
    ]
    colors = palette(:default)

    # Value of one method's column at one (dim, count): NaN if the file is missing
    function value_at(df, col)
        (df === nothing || !(col in propertynames(df))) && return NaN
        sigma === nothing && return mean(df[!, col])
        return df[argmin(abs.(df.Sigma .- sigma)), col]
    end

    subplots = Plots.Plot[]
    for (j, g) in enumerate(gaussians)
        dfs = Dict{Int, Any}()
        for d in dims
            csv_path = joinpath(input_dir, "$(name_prefix)_dim_$(d)_$(name_infix)_$(g).csv")
            if isfile(csv_path)
                dfs[d] = CSV.read(csv_path, DataFrame)
            else
                @warn "Missing results file" csv_path
                dfs[d] = nothing
            end
        end

        subtitle = sigma === nothing ? "mean over sigma" : "sigma = $sigma"
        p = plot(
            title = "$count_label: $g\n$subtitle",
            titlefontsize = 8,
            xlabel = "Dimension",
            ylabel = metric,
            guidefontsize = 7,
            tickfontsize = 6,
            xscale = :log2,
            xticks = (dims, string.(dims)),
            ylims = (0, 1),
            legend = false,
            framestyle = :box,
        )
        for (k, (prefix, label, marker)) in enumerate(methods)
            col = Symbol("$(prefix)_$(metric)")
            ys = [value_at(dfs[d], col) for d in dims]
            plot!(p, dims, ys; color = colors[k], marker = marker,
                  markersize = 3, linewidth = 1.5, label = label)
        end
        push!(subplots, p)
    end

    panels = plot(subplots...; layout = (1, length(gaussians)))

    # Legend-only panel shared by all subplots
    legend_plot = plot(; framestyle = :none, legend = :inside, legend_columns = -1,
                       legendfontsize = 7, xlims = (0, 1), ylims = (0, 1))
    for (k, (_, label, marker)) in enumerate(methods)
        plot!(legend_plot, [NaN], [NaN]; color = colors[k], marker = marker,
              markersize = 3, linewidth = 1.5, label = label)
    end

    fig = plot(panels, legend_plot;
               layout = grid(2, 1, heights = [0.86, 0.14]),
               size = (450 * length(gaussians), 480))

    mkpath(dirname(output_path))
    savefig(fig, output_path)
    println("Saved plot to $output_path")
    return fig
end

# Same plot for the ring data written by rings_sweep.jl.
function plot_rings_by_dimension(;
    dims::Vector{Int} = [4, 16, 64],
    ring_counts::Vector{Int} = [4, 8],
    metric::String = "NMI",
    sigma::Union{Real, Nothing} = nothing,
    input_dir::String = "data/Systhetic/",
    name_prefix::String = "final_rings_results",
)
    return plot_by_dimension(;
        dims = dims,
        gaussians = ring_counts,
        metric = metric,
        sigma = sigma,
        input_dir = input_dir,
        name_prefix = name_prefix,
        name_infix = "rings",
        count_label = "Rings",
        output_path = joinpath(input_dir, "rings_by_dim_$(metric).png"),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_by_dimension(metric = "NMI")
    plot_by_dimension(metric = "ARI")
    plot_rings_by_dimension(metric = "NMI")
    plot_rings_by_dimension(metric = "ARI")
end
