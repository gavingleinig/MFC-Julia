using MetricForestCompletion
using Random
using Clustering
using Plots

include("plot_utils_convex_hull.jl")

function run_and_plot_rings(; dim, num_rings, points_per_ring, radius_range, noise_sigma)
    PrecisionType = Float64
    num_clusters = num_rings
    dist_func = euclidean

    # --- Generate ring data instead of Gaussian blobs ---
    points, ground_truth = generate_rings(
        PrecisionType;
        dim = dim,
        num_rings = num_rings,
        points_per_ring = points_per_ring,
        radius_range = (PrecisionType(radius_range[1]), PrecisionType(radius_range[2])),
        noise_sigma = PrecisionType(noise_sigma),
    )

    # --- Same clustering calls as the Gaussian sweep / plotting script ---
    clustering_result = k_centering(points, num_clusters, 1, dist_func)
    points_matrix = stack(points)
    kmeans_result = kmeans(points_matrix, num_clusters)

    mfc_approx_result  = metric_forest_completion_approx(points, num_clusters, clustering_result.assignments, dist_func)
    mfc_optimal_result = metric_forest_completion_optimal(points, num_clusters, clustering_result.assignments, dist_func)
    mfc_simple_result  = metric_forest_completion_simple(points, num_clusters, clustering_result.assignments, dist_func)
    multi_rep_mfc_result = mulit_reps_metric_forest_completion(points, num_clusters, clustering_result.assignments, dist_func)

    all_ST_approx_edges  = vcat(mfc_approx_result.cluster_edges, mfc_approx_result.completion_edges)
    all_ST_optimal_edges = vcat(mfc_optimal_result.cluster_edges, mfc_optimal_result.completion_edges)
    all_ST_simple_edges  = vcat(mfc_simple_result.cluster_edges, mfc_simple_result.completion_edges)
    all_ST_multi_edges   = vcat(multi_rep_mfc_result.cluster_edges, multi_rep_mfc_result.completion_edges)

    n_pts = length(points)
    mst_complete      = convert_mst_to_weight(mst_implicit(points, dist_func))
    naive_st_complete = convert_mst_to_weight(naive_random_st(points, dist_func))

    best_mst, mst_ari         = best_single_linkage_threshold_increment(mst_complete, n_pts, ground_truth, .01)
    best_naive, naive_ari     = best_single_linkage_threshold_increment(naive_st_complete, n_pts, ground_truth, .01)
    best_approx, approx_ari   = best_single_linkage_threshold_increment(all_ST_approx_edges, n_pts, ground_truth, .01)
    best_optimal, optimal_ari = best_single_linkage_threshold_increment(all_ST_optimal_edges, n_pts, ground_truth, .01)
    best_simple, simple_ari   = best_single_linkage_threshold_increment(all_ST_simple_edges, n_pts, ground_truth, .01)
    best_multi, multi_ari     = best_single_linkage_threshold_increment(all_ST_multi_edges, n_pts, ground_truth, .01)

    kmeans_ari = randindex(kmeans_result.assignments, ground_truth)[1]
    kc_ari     = randindex(clustering_result.assignments, ground_truth)[1]

    # --- Project once, reuse for every panel ---
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
    outpath = joinpath(@__DIR__, "runings_clusters_vs_ground_truth_sigma_$(noise_sigma).png")
    savefig(combined, outpath)
    println("Saved: ", outpath)
    return combined
end

# ------------------------------------------------------------------
# Rings are a deliberately hard case for center-based methods: k-means and
# k-centering both assign points to the nearest *centroid*, but concentric
# rings share the same center -- so a centroid-based method has no way to
# tell an inner ring from an outer one by distance alone. Expect their ARI
# to be low here even though the rings are visually obvious. MST/graph-based
# methods that follow local neighbor structure should do much better, as
# long as noise_sigma doesn't make adjacent rings' points overlap.
# ------------------------------------------------------------------

#hi = 1 + (num_rings - 1) * 6 * noise_sigma 
run_and_plot_rings(
    dim = 2,
    num_rings = 16,
    points_per_ring = 300,
    radius_range = (2.0, 200.0),
    noise_sigma = .3,
)
