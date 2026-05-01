using Plots
using DataFrames
using CSV

function plot_clustering_results(df::DataFrame)

    sigmas = df.Sigma
    
    # ARI PLOT
    p_ari = plot(
        title = "Adjusted Rand Index (ARI) vs Sigma",
        xlabel = "Sigma",
        ylabel = "ARI",
        legend = :topright,
        grid = true
    )
    
    plot!(p_ari, sigmas, df.KC_ARI, label="K-Centering", marker=:circle, lw=2)
    plot!(p_ari, sigmas, df.MST_ARI, label="Optimal MST", marker=:square, lw=2)
    plot!(p_ari, sigmas, df.Naive_ARI, label="Naive ST", marker=:dtriangle, lw=2)
    plot!(p_ari, sigmas, df.MFC_Approx_ARI, label="MFC Approx", marker=:diamond, lw=2)
    plot!(p_ari, sigmas, df.MFC_Optimal_ARI, label="MFC Optimal", marker=:star5, lw=2)
    plot!(p_ari, sigmas, df.MFC_Simple_ARI, label="MFC Simple", marker=:cross, lw=2)
    plot!(p_ari, sigmas, df.KMeans_ARI, label="K-Means", marker=:xcross, lw=2)
    # NMI PLOT
    p_nmi = plot(
        title = "Normalized Mutual Info (NMI) vs Sigma",
        xlabel = "Sigma",
        ylabel = "NMI",
        legend = :topright,
        grid = true
    )
    
    plot!(p_nmi, sigmas, markersize=1, df.KC_NMI, label="K-Centering", marker=:circle, lw=2)
    plot!(p_nmi, sigmas, markersize=1, df.MST_NMI, label="Optimal MST", marker=:square, lw=2)
    plot!(p_nmi, sigmas, markersize=1, df.Naive_NMI, label="Naive ST", marker=:dtriangle, lw=2)
    plot!(p_nmi, sigmas, markersize=1, df.MFC_Approx_NMI, label="MFC Approx", marker=:diamond, lw=2)
    plot!(p_nmi, sigmas, markersize=1, df.MFC_Optimal_NMI, label="MFC Optimal", marker=:star5, lw=2)
    plot!(p_nmi, sigmas, markersize=1, df.MFC_Simple_NMI, label="MFC Simple", marker=:cross, lw=2)
    plot!(p_nmi, sigmas, markersize=1, df.KMeans_NMI, label="K-Means", marker=:xcross, lw=2)


    final_plot = plot(p_ari, p_nmi, layout=(1, 2), size=(1000, 450), margin=5Plots.mm)
    
    # Display the plot
    display(final_plot)
    
    output_image = "clustering_metrics_vs_sigma.png"
    savefig(final_plot, output_image)
    println("Plot saved to $output_image")
    
    return final_plot
end

df= CSV.read("data/Systhetic/clustering_results.csv", DataFrame)
plot_clustering_results(df)