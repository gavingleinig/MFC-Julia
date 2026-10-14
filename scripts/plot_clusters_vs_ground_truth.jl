using MetricForestCompletion
using Random
using Clustering
using Plots

include("plot_utils_convex_hull.jl")

function run_and_plot_single_sigma(; dim, num_gauss, num_points, mean_range, sigma)
    PrecisionType = Float64
    num_clusters = num_gauss
    dist_func = euclidean

    # --- Same generation call as run_gaussian_sigma_sweep, one fixed sigma ---
    points, ground_truth = generate_gaussians(
        PrecisionType;
        dim = dim,
        num_centers = num_gauss,
        points_per_center = num_points ÷ num_gauss,
        mean_range = (PrecisionType(mean_range[1]), PrecisionType(mean_range[2])),
        sigma_range = (PrecisionType(sigma), PrecisionType(sigma)),
    )

    # --- Same clustering calls as the sweep ---
    clustering_result = k_centering(points, num_clusters, 1, dist_func)
    points_matrix = stack(points)
    #kmeans_result = kmeans(points_matrix, num_clusters)

    mfc_approx_result  = metric_forest_completion_approx(points, num_clusters, clustering_result.assignments, dist_func)
    mfc_optimal_result = metric_forest_completion_optimal(points, num_clusters, clustering_result.assignments, dist_func)
    mfc_simple_result  = metric_forest_completion_simple(points, num_clusters, clustering_result.assignments, dist_func)
    multi_rep_mfc_result = mulit_reps_metric_forest_completion(points, num_clusters, clustering_result.assignments, dist_func)

    all_ST_approx_edges  = vcat(mfc_approx_result.cluster_edges, mfc_approx_result.completion_edges)
    all_ST_optimal_edges = vcat(mfc_optimal_result.cluster_edges, mfc_optimal_result.completion_edges)
    all_ST_simple_edges  = vcat(mfc_simple_result.cluster_edges, mfc_simple_result.completion_edges)
    all_ST_multi_edges   = vcat(multi_rep_mfc_result.cluster_edges, multi_rep_mfc_result.completion_edges)

    n_pts = length(points)
    mst_complete       = convert_mst_to_weight(mst_implicit(points, dist_func))
    #naive_st_complete  = convert_mst_to_weight(naive_random_st(points, dist_func))

    best_mst, mst_ari         = best_single_linkage_threshold_increment(mst_complete, n_pts, ground_truth, .01)
    #best_naive, naive_ari     = best_single_linkage_threshold_increment(naive_st_complete, n_pts, ground_truth, .01)
    best_approx, approx_ari   = best_single_linkage_threshold_increment(all_ST_approx_edges, n_pts, ground_truth, .01)
    best_optimal, optimal_ari = best_single_linkage_threshold_increment(all_ST_optimal_edges, n_pts, ground_truth, .01)
    best_simple, simple_ari   = best_single_linkage_threshold_increment(all_ST_simple_edges, n_pts, ground_truth, .01)
    best_multi, multi_ari     = best_single_linkage_threshold_increment(all_ST_multi_edges, n_pts, ground_truth, .01)

    # KMeans / K-Centering don't come with an ARI from the single-linkage
    # search above (that's only run on tree/edge-based methods), so compute
    # it the same way the sweep script does: randindex(...)[1] is the ARI.
    #kmeans_ari = randindex(kmeans_result.assignments, ground_truth)[1]
    kc_ari     = randindex(clustering_result.assignments, ground_truth)[1]

    # --- Project once, reuse the same 2D coords for every panel ---
    coords2d = to_2d(points)

    methods = [
        #("KMeans",       kmeans_result.assignments,      kmeans_ari),
        ("K-Centering",  clustering_result.assignments,  kc_ari),
        ("MST",          best_mst.assignments,           mst_ari),
        #("Naive ST",     best_naive.assignments,         naive_ari),
        ("MFC Approx",   best_approx.assignments,        approx_ari),
        ("MFC Optimal",  best_optimal.assignments,       optimal_ari),
        ("MFC Simple",   best_simple.assignments,        simple_ari),
        ("MFC Multi",    best_multi.assignments,         multi_ari),
    ]

    panels = [
        plot_clusters_vs_ground_truth(
            coords2d, ground_truth, assignments;
            title = "$name (ARI = $(round(ari, digits=3)))",
        )
        for (name, assignments, ari) in methods
    ]

    combined = plot(panels...; layout = (4, 2), size = (1300, 2000))
    outpath = joinpath(@__DIR__, "clusters_vs_ground_truth_sigma_$(sigma).png")
    savefig(combined, outpath)
    println("Saved: ", outpath)
    return combined
end

# ------------------------------------------------------------------
# Example usage.
#
# NOTE: the real sweep uses dim=16, num_gauss=128, num_points=2000 — that's
# 128 hulls per panel, which will look extremely dense. Start with something
# small to sanity-check the plotting itself, then scale up to match the
# real experiment.
# ------------------------------------------------------------------
run_and_plot_single_sigma(
    dim = 2,
    num_gauss = 6,
    num_points = 600,
    mean_range = (-5, 5),
    sigma = 1.5,
)

# To reproduce a panel from the actual sweep, e.g.:
# run_and_plot_single_sigma(dim=16, num_gauss=128, num_points=2000, mean_range=(-5,5), sigma=0.5)