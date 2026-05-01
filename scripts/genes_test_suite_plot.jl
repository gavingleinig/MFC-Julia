using Plots

taxonomies = [
    (2, "Phylum", 3),
    (3, "Class", 20),
    (4, "Order", 38),
    (5, "Family", 67),
    (6, "Genus", 130),
    (7, "Species", 148)
]

x_values = [t[3] for t in taxonomies]


x_labels = ["$(t[2])\n($(t[3]))" for t in taxonomies] 
my_xticks = (x_values, x_labels)

ari_init_k     = [0.1722, 0.7225, 0.7145, 0.8015, 0.7895, 0.7607]
ari_mst        = [0.3475, 0.9694, 0.7517, 0.7935, 0.8043, 0.8163]
ari_naive      = [0.0055, 0.0083, 0.0063, 0.0029, 0.0056, 0.0055]
ari_mfc_approx = [0.3213, 0.9552, 0.8268, 0.8541, 0.8119, 0.8135]
ari_mfc_opt    = [0.3078, 0.9359, 0.8232, 0.8267, 0.8040, 0.8159]
ari_mfc_simp   = [0.2869, 0.8657, 0.8148, 0.8545, 0.8281, 0.8013]

nmi_init_k     = [0.1645, 0.7763, 0.8167, 0.8711, 0.8689, 0.8670]
nmi_mst        = [0.5041, 0.9239, 0.8504, 0.8543, 0.8585, 0.8686]
nmi_naive      = [0.0324, 0.0936, 0.1184, 0.3171, 0.2469, 0.2624]
nmi_mfc_approx = [0.4733, 0.9059, 0.8805, 0.8922, 0.8901, 0.8661]
nmi_mfc_opt    = [0.4631, 0.9090, 0.8751, 0.8772, 0.8577, 0.8679]
nmi_mfc_simp   = [0.4534, 0.8659, 0.8666, 0.8897, 0.8771, 0.8702]

l = @layout [grid(1, 2)
             b{0.1h}] 

final_plot = plot(
    layout = l,
    size = (1100, 650),
    bottom_margin = 5Plots.mm,
    left_margin = 10Plots.mm,
    grid = true
)

plot!(final_plot, subplot = 1,
    title = "ARI Score vs Taxonomy",
    ylabel = "ARI Score",
    legend = false,
    xticks = my_xticks,
    xrotation = 45,
    xflip = true,
    lw = 1
)

plot!(final_plot, x_values, ari_init_k,     subplot=1, marker=:circle,    lw=1, markersize=3)
plot!(final_plot, x_values, ari_mst,        subplot=1, marker=:square,    lw=1, markersize=3)
plot!(final_plot, x_values, ari_naive,      subplot=1, marker=:dtriangle, lw=1, markersize=3)
plot!(final_plot, x_values, ari_mfc_approx, subplot=1, marker=:diamond,   lw=1, markersize=3)
plot!(final_plot, x_values, ari_mfc_opt,    subplot=1, marker=:star5,     lw=1, markersize=3)
plot!(final_plot, x_values, ari_mfc_simp,   subplot=1, marker=:cross,     lw=1, markersize=3)

plot!(final_plot, subplot = 2,
    title = "NMI Score vs Taxonomy",
    ylabel = "NMI Score",
    legend = false,
    xticks = my_xticks,
    xrotation = 45,
    xflip = true,
    lw = 1
)

plot!(final_plot, x_values, nmi_init_k,     subplot=2, marker=:circle,    lw=1, markersize=3)
plot!(final_plot, x_values, nmi_mst,        subplot=2, marker=:square,    lw=1, markersize=3)
plot!(final_plot, x_values, nmi_naive,      subplot=2, marker=:dtriangle, lw=1, markersize=3)
plot!(final_plot, x_values, nmi_mfc_approx, subplot=2, marker=:diamond,   lw=1, markersize=3)
plot!(final_plot, x_values, nmi_mfc_opt,    subplot=2, marker=:star5,     lw=1, markersize=3)
plot!(final_plot, x_values, nmi_mfc_simp,   subplot=2, marker=:cross,     lw=1, markersize=3)


plot!(final_plot, rand(0), subplot = 3, 
    label = "Init K-Centering", 
    marker = :circle,    
    legend = :inside, 
    legend_column = 6, 
    foreground_color_legend = :transparent,
    background_color_legend = :transparent, 
    framestyle = :none, 
    axis = false, 
    grid = false
)
plot!(final_plot, rand(0), subplot = 3, label = "MST Single",       marker = :square,    legend = :inside)
plot!(final_plot, rand(0), subplot = 3, label = "Naive ST",         marker = :dtriangle, legend = :inside)
plot!(final_plot, rand(0), subplot = 3, label = "MFC Approx",       marker = :diamond,   legend = :inside)
plot!(final_plot, rand(0), subplot = 3, label = "MFC Optimal",      marker = :star5,     legend = :inside)
plot!(final_plot, rand(0), subplot = 3, label = "MFC Simple",       marker = :cross,     legend = :inside)

display(final_plot)

savefig(final_plot, "taxonomy_real_scores_flipped_matched.png")