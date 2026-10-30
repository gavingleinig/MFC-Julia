using CSV
using DataFrames
using Plots

# Plots every method's ARI (or NMI) against sigma on a single graph,
# from a CSV written by run_gaussian_sigma_sweep in evaluation.jl.
function plot_sweep(csv_path::String;
    metric::String = "ARI",             # "ARI" or "NMI"
    output_path::Union{String, Nothing} = nothing,
    title::String = basename(csv_path)
)
    df = CSV.read(csv_path, DataFrame)

    # (column prefix, legend label, is MFC variant)
    methods = [
        ("KMeans",      "KMeans",      false),
        ("KC",          "K-Centering", false),
        ("MST",         "MST",         false),
        ("Naive",       "Naive ST",    false),
        ("MFC_Approx",  "MFC Approx",  true),
        ("MFC_Optimal", "MFC Optimal", true),
        ("MFC_Simple",  "MFC Simple",  true),
        ("MFC_Multi",   "MFC Multi",   true),
    ]

    plt = plot(
        xlabel = "Sigma",
        ylabel = metric,
        title = title,
        legend = :outertopright,
        ylims = (0, 1.05),
        size = (950, 520),
    )

    for (prefix, label, is_mfc) in methods
        col = Symbol("$(prefix)_$(metric)")
        col in propertynames(df) || continue   # skip methods missing from the CSV
        plot!(plt, df.Sigma, df[!, col];
              label = label,
              linewidth = 2,
              marker = :circle,
              markersize = 4,
              linestyle = is_mfc ? :dash : :solid)
    end

    if output_path === nothing
        output_path = replace(csv_path, r"\.csv$" => "") * "_$(metric).png"
    end
    savefig(plt, output_path)
    println("Saved plot to $output_path")

    return plt
end

if abspath(PROGRAM_FILE) == @__FILE__
    csv_path = isempty(ARGS) ?
        "data/Systhetic/final_clustering_results_dim_16_gaus_128.csv" : ARGS[1]
    plot_sweep(csv_path; metric = "ARI")
    plot_sweep(csv_path; metric = "NMI")
end
