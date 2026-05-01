using MetricForestCompletion
using Random
using Clustering
using DataFrames
using CSV

function run_gaussian_sigma_sweep(sigma_values::AbstractVector,in_dim::Int64,in_gauss::Int64,CSV_label::String )
    # Setup Params
    PrecisionType = Float64
    
    dim = in_dim
    num_gauss = in_gauss
    num_points = 20000
    num_clusters = in_gauss 
    
    # Distance Metric
    dist_func = euclidean 

    results = []

    results_runtime = []

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
        k_centering_runtime = @elapsed clustering_result = k_centering(points, num_clusters, 1, dist_func)

        points_matrix = stack(points) # put points into 2D vector

        
        full_kmeans_runtime = @elapsed kmeans_result = kmeans(points_matrix, num_clusters)
        

        # MFC Variants
        mfc_approx_runtime  = @elapsed mfc_approx_result = metric_forest_completion_approx(points, num_clusters, clustering_result.assignments, dist_func)
        mfc_optimal_runtime = @elapsed mfc_optimal_result = metric_forest_completion_optimal(points, num_clusters, clustering_result.assignments, dist_func)
        mfc_simple_runtime  = @elapsed mfc_simple_result = metric_forest_completion_simple(points, num_clusters, clustering_result.assignments, dist_func)

        # Combine all edges 
        all_ST_approx_edges  = vcat(mfc_approx_result.cluster_edges, mfc_approx_result.completion_edges)
        all_ST_optimal_edges = vcat(mfc_optimal_result.cluster_edges, mfc_optimal_result.completion_edges)
        all_ST_simple_edges  = vcat(mfc_simple_result.cluster_edges, mfc_simple_result.completion_edges)
        
        n_pts = length(points) 

        # Compute Optimal MST & Naive ST
        mst_runtime      = @elapsed mst_complete  = convert_mst_to_weight(mst_implicit(points, dist_func))
        naive_st_runtime = @elapsed naive_st_complete = convert_mst_to_weight(naive_random_st(points, dist_func))

        
        gamma_overlap = compute_gamma_overlap(
            mfc_approx_result.cluster_edges, 
            mst_complete, 
            clustering_result.assignments
        )

        mst_single_link_runtime     = @elapsed best_mst, mst_ari         = best_single_linkage_threshold_increment(mst_complete, n_pts, ground_truth,.01 )
        naive_single_link_runtime   = @elapsed best_naive, naive_ari     = best_single_linkage_threshold_increment(naive_st_complete, n_pts, ground_truth,.01 )
        approx_single_link_runtime  = @elapsed best_approx, approx_ari   = best_single_linkage_threshold_increment(all_ST_approx_edges, n_pts, ground_truth, .01 )
        optimal_single_link_runtime = @elapsed best_optimal, optimal_ari = best_single_linkage_threshold_increment(all_ST_optimal_edges, n_pts, ground_truth, .01 )
        simple_single_link_runtime  = @elapsed best_simple, simple_ari   = best_single_linkage_threshold_increment(all_ST_simple_edges, n_pts, ground_truth, .01 )


        # Full Runtime Caluculations
        full_mst_runtime     = mst_runtime                                 #  + mst_single_link_runtime
        full_naive_runtime   = naive_st_runtime                            #  + naive_single_link_runtime 
        full_approx_runtime  = mfc_approx_runtime + k_centering_runtime    # + approx_single_link_runtime  
        full_optimal_runtime = mfc_optimal_runtime+ k_centering_runtime    # + optimal_single_link_runtime + k_centering_runtime
        full_simple_runtime  = mfc_simple_runtime + k_centering_runtime    # + simple_single_link_runtime  + k_centering_runtime
        println("mfc_approx_runtime cluster runtime: ",mfc_approx_result.sub_cluster_runtime)
        println("mfc_approx_runtime  completion runtime: ",mfc_approx_result.completion_edges_runtime )
        println("mfc_approx_runtime  full runtime: ",mfc_approx_result.completion_runtime )
        println("mfc_optimal_runtime cluster runtime: ",mfc_optimal_result.sub_cluster_runtime)
        println("mfc_optimal_runtime  completion runtime: ",mfc_optimal_result.completion_edges_runtime )
        println("mfc_optimal_runtime  full runtime: ",mfc_optimal_result.completion_runtime )





        # Clustering results
        push!(results, (
            Sigma           = sigma,

            Gamma_Overlap   = gamma_overlap,

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

        # Runtime results
        push!(results_runtime, (
            Sigma               = sigma,

            KMeans_runtime      = full_kmeans_runtime,

            KC_runtime          = k_centering_runtime,
            
            MST_runtime         = full_mst_runtime,
            
            Naive_runtime       = full_naive_runtime,
            
            MFC_Approx_runtime  = full_approx_runtime,
            
            MFC_Optimal_runtime = full_optimal_runtime,
            
            MFC_Simple_runtime  = full_simple_runtime,
        ))



        loop_elapsed = time() - loop_start
        println("   Completed in $(round(loop_elapsed, digits=2)) sec")
    end

    # convert results to a DataFrame, print
    df = DataFrame(results)

    runtime_df =  DataFrame(results_runtime)
    
    println("\nFinal Results:")
    println(df)

    println("\nRuntime Results:")
    println(runtime_df)

    # Save to CSV
    output_file = CSV_label
    runtime_output_file = string("runtime_",CSV_label)


    CSV.write(output_file, df)
    println("\nResults saved to $output_file")

    CSV.write(runtime_output_file, runtime_df)
    println("\nRuntime safed to $runtime_output_file")
    
    return df
end

sigmas_to_test = 0.1:0.2:2
df_results = run_gaussian_sigma_sweep(sigmas_to_test,256,128,"clustering_results.csv")