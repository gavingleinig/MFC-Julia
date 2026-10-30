using CSV
using DataFrames
using Plots

const METHOD_LABELS = Dict(
    "KMeans"      => "KMeans",
    "KC"          => "K-Centering",
    "MST"         => "Optimal MST",
    "Naive"       => "Naive ST",
    "MFC_Approx"  => "MFC Approx",
    "MFC_Optimal" => "MFC Optimal",
    "MFC_Simple"  => "MFC Simple",
    "MFC_Multi"   => "MFC Multi",
)

# One panel: mean ARI and mean NMI per method as side-by-side bars,
# with +/- one standard deviation (over the trials) as error bars.
function summary_panel(csv_path::String, panel_title::String)
    df = CSV.read(csv_path, DataFrame)
    labels = [get(METHOD_LABELS, m, m) for m in df.Method]
    x = collect(1:nrow(df))

    p = plot(
        title = panel_title,
        titlefontsize = 9,
        ylabel = "Score (mean ± std over trials)",
        guidefontsize = 7,
        tickfontsize = 7,
        xticks = (x, labels),
        xrotation = 35,
        ylims = (0, 1.1),
        legend = :topright,
        legendfontsize = 7,
        framestyle = :box,
    )
    bar!(p, x .- 0.2, df.Mean_ARI; yerror = df.Std_ARI, bar_width = 0.35,
         label = "ARI", color = palette(:default)[1])
    bar!(p, x .+ 0.2, df.Mean_NMI; yerror = df.Std_NMI, bar_width = 0.35,
         label = "NMI", color = palette(:default)[2])
    return p
end

# Plots the trial summaries (from gaussian_trials.jl / rings_trials.jl)
# side by side in one figure. Panels whose file is missing are skipped.
function plot_summary(;
    input_dir::String = "data/Systhetic/",
    gaussian_csv::String = joinpath(input_dir, "gaussian_trials_summary.csv"),
    rings_csv::String = joinpath(input_dir, "rings_trials_summary.csv"),
    output_path::String = joinpath(input_dir, "trials_summary.png"),
)
    panels = Plots.Plot[]
    for (path, name) in ((gaussian_csv, "Gaussians"), (rings_csv, "Rings"))
        if isfile(path)
            push!(panels, summary_panel(path, name))
        else
            @warn "Missing summary file, skipping" path
        end
    end
    isempty(panels) && error("No summary files found")

    fig = plot(panels...; layout = (1, length(panels)),
               size = (520 * length(panels), 480), margin = 4Plots.mm)

    mkpath(dirname(output_path))
    savefig(fig, output_path)
    println("Saved plot to $output_path")
    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_summary()
end
