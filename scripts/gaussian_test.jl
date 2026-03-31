using MetricForestCompletion
using Random

function run_gaussian_test()
    # Setup Params
    PrecisionType = Float64
    
    dim = 4
    num_gauss = 100
    points_per_gauss = 200
    clusters = 16 
    
    # Distance Metric
    dist_func = euclidean 

    # Generate Dataset
    points = generate_gaussians(
        PrecisionType; 
        dim = dim, 
        num_centers = num_gauss, 
        points_per_center = points_per_gauss,
        mean_range = (PrecisionType(-5.0), PrecisionType(5.0)),
        sigma_range = (PrecisionType(0.5), PrecisionType(0.8))
    )

    # Initial Clustering
    clustering_result = k_centering(points, clusters, dist_func)

    # MFC 
    mfc_result = metric_forest_completion(
        points, 
        clusters, 
        clustering_result.assignments,
        dist_func
    )

    # Print results
    println("Test Results")
    println("Data Type:           ", eltype(eltype(points)))
    println("Points Generated:    ", length(points))
    println("Clustering Runtime:  ", clustering_result.runtime, " seconds")
    println("MFC Total Runtime:   ", mfc_result.completion_runtime, " seconds")
    println("Cluster Edges Found: ", length(mfc_result.cluster_edges))
end

run_gaussian_test()