using CSV
using DataFrames
using Plots

# (summary "Method" value, legend label, marker, dashed?)
const SUMMARY_METHODS = [
    ("KMeans",      "KMeans",      :utriangle, false),
    ("KC",          "K-Centering", :circle,    false),
    ("MST",         "Optimal MST", :square,    false),
    ("Naive",       "Naive ST",    :dtriangle, false),
    ("MFC_Approx",  "MFC Approx",  :diamond,   true),
    ("MFC_Optimal", "MFC Optimal", :hexagon,   true),
    ("MFC_Simple",  "MFC Simple",  :xcross,    true),
    ("MFC_Multi",   "MFC Multi",   :star5,     true),
]

# Finds files named <name_prefix>_dim_<d>_gaus_<g>_summary(.csv) in `input_dir`.
# Returns Dict((dim, gaussians) => path).
function find_summary_files(input_dir::String, name_prefix::String)
    pattern = Regex("^" * name_prefix * "_dim_(\\d+)_gaus_(\\d+)_summary(\\.csv)?\$")
    found = Dict{Tuple{Int, Int}, String}()
    for f in readdir(input_dir)
        m = match(pattern, f)
        m === nothing && continue
        found[(parse(Int, m.captures[1]), parse(Int, m.captures[2]))] = joinpath(input_dir, f)
    end
    return found
end

# Accuracy vs. dimension from the per-configuration summary tables
# (Method, Mean_ARI, Std_ARI, Mean_NMI, Std_NMI, Mean_runtime).
# One column of panels per Gaussian count, one row per metric, every method as a
# line with +/- std error bars, and one shared legend underneath.
# `dims` / `gaussians` default to whatever is found in `input_dir`.
function plot_summary_by_dimension(;
    input_dir::String = "data/Systhetic/",
    name_prefix::String = "final_clustering_results",
    metrics::Vector{String} = ["ARI", "NMI"],
    dims::Union{Vector{Int}, Nothing} = nothing,
    gaussians::Union{Vector{Int}, Nothing} = nothing,
    output_path::String = joinpath(input_dir, "summary_by_dim.png"),
)
    files = find_summary_files(input_dir, name_prefix)
    isempty(files) && error("No files like $(name_prefix)_dim_<d>_gaus_<g>_summary found in $input_dir")

    dims === nothing && (dims = sort(unique(first.(keys(files)))))
    gaussians === nothing && (gaussians = sort(unique(last.(keys(files)))))

    colors = palette(:default)

    # Mean and std of one method/metric at one (dim, gaussians); NaN if missing
    function lookup(df, method, metric)
        df === nothing && return (NaN, 0.0)
        row = df[df.Method .== method, :]
        isempty(row) && return (NaN, 0.0)
        return (row[1, Symbol("Mean_$metric")], row[1, Symbol("Std_$metric")])
    end

    subplots = Plots.Plot[]
    for (i, metric) in enumerate(metrics), (j, g) in enumerate(gaussians)
        p = plot(
            title = i == 1 ? "Gaussians: $g" : "",
            titlefontsize = 9,
            xlabel = i == length(metrics) ? "Dimension" : "",
            ylabel = j == 1 ? metric : "",
            guidefontsize = 8,
            tickfontsize = 7,
            xscale = :log2,
            xticks = (dims, string.(dims)),
            ylims = (0, 1.05),
            legend = false,
            framestyle = :box,
        )

        dfs = Dict(d => (haskey(files, (d, g)) ? CSV.read(files[(d, g)], DataFrame) : nothing)
                   for d in dims)
        any(isnothing, values(dfs)) && @warn "Some summary files are missing for gaussians = $g"

        for (k, (method, label, marker, dashed)) in enumerate(SUMMARY_METHODS)
            vals = [lookup(dfs[d], method, metric) for d in dims]
            plot!(p, dims, first.(vals); yerror = last.(vals),
                  color = colors[k], marker = marker, markersize = 4, linewidth = 1.5,
                  linestyle = dashed ? :dash : :solid, label = label)
        end
        push!(subplots, p)
    end

    panels = plot(subplots...; layout = (length(metrics), length(gaussians)))

    # Legend-only panel shared by all subplots
    legend_plot = plot(; framestyle = :none, legend = :inside, legend_columns = -1,
                       legendfontsize = 8, xlims = (0, 1), ylims = (0, 1))
    for (k, (_, label, marker, dashed)) in enumerate(SUMMARY_METHODS)
        plot!(legend_plot, [NaN], [NaN]; color = colors[k], marker = marker, markersize = 4,
              linewidth = 1.5, linestyle = dashed ? :dash : :solid, label = label)
    end

    fig = plot(panels, legend_plot;
               layout = grid(2, 1, heights = [0.9, 0.1]),
               size = (480 * length(gaussians), 340 * length(metrics) + 70))

    mkpath(dirname(output_path))
    savefig(fig, output_path)
    println("Saved plot to $output_path")
    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_summary_by_dimension()
end
