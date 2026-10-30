using CSV
using DataFrames
using Plots

# Grid of sigma-sweep plots: one row per dimension, one column per Gaussian count,
# a single shared legend underneath. Reads the CSVs written by dimension_sweep.jl
# (named <name_prefix>_dim_<d>_gaus_<g>.csv).
function plot_sweep_grid(;
    dims::Vector{Int} = [4, 16, 64],
    gaussians::Vector{Int} = [16, 64],
    metric::String = "NMI",                     # "ARI" or "NMI"
    input_dir::String = "data/Systhetic/",
    name_prefix::String = "final_clustering_results",
    output_path::String = joinpath(input_dir, "sweep_grid_$(metric).png"),
)
    # (CSV column prefix, legend label, marker). Add ("KMeans", "KMeans", :utriangle)
    # or ("MFC_Multi", "MFC Multi", :star5) here to include them.
    methods = [
        ("KC",          "K-Centering", :circle),
        ("MST",         "Optimal MST", :square),
        ("Naive",       "Naive ST",    :dtriangle),
        ("MFC_Approx",  "MFC Approx",  :diamond),
        ("MFC_Optimal", "MFC Optimal", :hexagon),
        ("MFC_Simple",  "MFC Simple",  :xcross),
    ]
    colors = palette(:default)

    subplots = Plots.Plot[]
    for (i, d) in enumerate(dims), (j, g) in enumerate(gaussians)
        p = plot(
            title = "Gaussians: $g\nd=$d",
            titlefontsize = 8,
            xlabel = i == length(dims) ? "Sigma" : "",
            ylabel = metric,
            guidefontsize = 7,
            tickfontsize = 6,
            ylims = (0, 1),
            legend = false,
            framestyle = :box,
        )

        csv_path = joinpath(input_dir, "$(name_prefix)_dim_$(d)_gaus_$(g).csv")
        if isfile(csv_path)
            df = CSV.read(csv_path, DataFrame)
            for (k, (prefix, label, marker)) in enumerate(methods)
                col = Symbol("$(prefix)_$(metric)")
                col in propertynames(df) || continue
                plot!(p, df.Sigma, df[!, col];
                      color = colors[k], marker = marker, markersize = 3,
                      linewidth = 1.5, label = label)
            end
        else
            @warn "Missing results file" csv_path
            annotate!(p, 0.5, 0.5, text("missing", 7))
        end
        push!(subplots, p)
    end

    grid_plot = plot(subplots...; layout = (length(dims), length(gaussians)))

    # Legend-only panel shared by all subplots
    legend_plot = plot(; framestyle = :none, legend = :inside, legend_columns = -1,
                       legendfontsize = 7, xlims = (0, 1), ylims = (0, 1))
    for (k, (_, label, marker)) in enumerate(methods)
        plot!(legend_plot, [NaN], [NaN]; color = colors[k], marker = marker,
              markersize = 3, linewidth = 1.5, label = label)
    end

    fig = plot(grid_plot, legend_plot;
               layout = grid(2, 1, heights = [0.92, 0.08]),
               size = (1000, 800))

    mkpath(dirname(output_path))
    savefig(fig, output_path)
    println("Saved plot to $output_path")
    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_sweep_grid(metric = "NMI")
    plot_sweep_grid(metric = "ARI")
end
