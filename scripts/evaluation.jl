using MetricForestCompletion
using Random
using Clustering
using DataFrames
using CSV

function run_gaussian_sigma_sweep(sigma_values::AbstractVector)
    # Setup Params
    PrecisionType = Float64
    
    dim = 8
    num_gauss = 32
    num_points = 20000
    num_clusters = 16 
    
    # Distance Metric
    dist_func = euclidean 

    results = []

    for sigma in sigma_values
        println("Evaluating sigma = $sigma")
        loop_start = time()

        # Generate Dataset
        points, ground_truth = generate_gaussians(
            PrecisionType; 
            dim = dim, 
            num_centers = num_gauss, 
            points_per_center = num_points ÷ num_gauss,
            mean_range = (PrecisionType(-5), PrecisionType(5)),
            sigma_range = (PrecisionType(sigma), PrecisionType(sigma))
        )

        # Initial Clustering
        clustering_result = k_centering(points, num_clusters, 1, dist_func)

        points_matrix = stack(points) # put points into 2D vector
        kmeans_result = kmeans(points_matrix, num_clusters)

        # MFC Variants
        mfc_approx_result = metric_forest_completion_approx(points, num_clusters, clustering_result.assignments, dist_func)
        mfc_optimal_result = metric_forest_completion_optimal(points, num_clusters, clustering_result.assignments, dist_func)
        mfc_simple_result = metric_forest_completion_simple(points, num_clusters, clustering_result.assignments, dist_func)

        # Combine all edges 
        all_ST_approx_edges  = vcat(mfc_approx_result.cluster_edges, mfc_approx_result.completion_edges)
        all_ST_optimal_edges = vcat(mfc_optimal_result.cluster_edges, mfc_optimal_result.completion_edges)
        all_ST_simple_edges  = vcat(mfc_simple_result.cluster_edges, mfc_simple_result.completion_edges)
        
        n_pts = length(points) 

        # Compute Optimal MST & Naive ST
        mst_complete = convert_mst_to_weight(mst_implicit(points, dist_func))
        naive_st_complete = convert_mst_to_weight(naive_random_st(points, dist_func))

        best_mst, mst_ari         = best_single_linkage_threshold_increment(mst_complete, n_pts, ground_truth,0.01 )
        best_naive, naive_ari     = best_single_linkage_threshold_increment(naive_st_complete, n_pts, ground_truth,0.01 )
        best_approx, approx_ari   = best_single_linkage_threshold_increment(all_ST_approx_edges, n_pts, ground_truth,0.01 )
        best_optimal, optimal_ari = best_single_linkage_threshold_increment(all_ST_optimal_edges, n_pts, ground_truth,0.01 )
        best_simple, simple_ari   = best_single_linkage_threshold_increment(all_ST_simple_edges, n_pts, ground_truth ,0.01)

        push!(results, (
            Sigma           = sigma,

            KMeans_ARI      = randindex(kmeans_result.assignments, ground_truth)[1],
            KMeans_NMI      = mutualinfo(kmeans_result.assignments, ground_truth, normed=true),

            KC_ARI          = randindex(clustering_result.assignments, ground_truth)[1],
            KC_NMI          = mutualinfo(clustering_result.assignments, ground_truth, normed=true),
            
            MST_ARI         = mst_ari,
            MST_NMI         = mutualinfo(best_mst.assignments, ground_truth, normed=true),
            
            Naive_ARI       = naive_ari,
            Naive_NMI       = mutualinfo(best_naive.assignments, ground_truth, normed=true),
            
            MFC_Approx_ARI  = approx_ari,
            MFC_Approx_NMI  = mutualinfo(best_approx.assignments, ground_truth, normed=true),
            
            MFC_Optimal_ARI = optimal_ari,
            MFC_Optimal_NMI = mutualinfo(best_optimal.assignments, ground_truth, normed=true),
            
            MFC_Simple_ARI  = simple_ari,
            MFC_Simple_NMI  = mutualinfo(best_simple.assignments, ground_truth, normed=true)
        ))
        loop_elapsed = time() - loop_start
        println("   Completed in $(round(loop_elapsed, digits=2)) sec")
    end

    # convert results to a DataFrame, print
    df = DataFrame(results)
    
    println("\nFinal Results:")
    println(df)

    # Save to CSV
    output_file = "clustering_results.csv"

    CSV.write(output_file, df)
    println("\nResults saved to $output_file")
    
    return df
end

sigmas_to_test = 0.1:0.2:1.5
df_results = run_gaussian_sigma_sweep(sigmas_to_test)