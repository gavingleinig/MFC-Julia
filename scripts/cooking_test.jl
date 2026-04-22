using MetricForestCompletion
using Random
using Clustering

function run_cooking_test(file_path::String)

    num_clusters = 16
    dist_func = jaccard

    points = Set{Int}[]
    edge_size_filter = 0 # Potentially filter out small sets
    
    # Load dataset from txt file; one set per line of comma seperated integers
    for line in eachline(file_path)
        parts = split(strip(line), ',')

        if isempty(parts) || parts == [""]
            continue
        end
        
        s = Set{Int}(parse(Int, p) for p in parts)
        
        if length(s) >= edge_size_filter
            push!(points, s)
        end
    end
    
    num_points = length(points)
    println("Loaded $(num_points) sets from cooking dataset.")

    # Initial Clustering
    clustering_result = k_centering(points, num_clusters, 1, dist_func)

    # MFC Variants
    mfc_approx_result = metric_forest_completion_approx(
        points, num_clusters, clustering_result.assignments, dist_func
    )

    mfc_optimal_result = metric_forest_completion_optimal(
        points, num_clusters, clustering_result.assignments, dist_func
    )

    mfc_simple_result = metric_forest_completion_simple(
        points, num_clusters, clustering_result.assignments, dist_func
    )

    all_ST_approx_edges  = vcat(mfc_approx_result.cluster_edges, mfc_approx_result.completion_edges)
    all_ST_optimal_edges = vcat(mfc_optimal_result.cluster_edges, mfc_optimal_result.completion_edges)
    all_ST_simple_edges  = vcat(mfc_simple_result.cluster_edges, mfc_simple_result.completion_edges)
    
    t_mst_compute = @elapsed mst_edges = mst_implicit(points, dist_func)
    mst_complete = convert_mst_to_weight(mst_edges)

    naive_st_edges = naive_random_st(points, dist_func)
    naive_st_complete = convert_mst_to_weight(naive_st_edges)



    println("Optimal MST Computation:             ", round(t_mst_compute, digits=4), " seconds")

    println("\n### MFC Runtimes ###")
    
    println("Initial Clustering:                  ", round(clustering_result.runtime, digits=4), " seconds")
    
    println("\nMFC Approx Sub-cluster Runtime:      ", round(mfc_approx_result.sub_cluster_runtime, digits=4), " seconds")
    println("MFC Approx Completion Edges Runtime: ", round(mfc_approx_result.completion_edges_runtime, digits=4), " seconds")
    println("MFC Approx Total Runtime:            ", round(mfc_approx_result.completion_runtime, digits=4), " seconds")

    println("\nMFC Optimal Sub-cluster Runtime:     ", round(mfc_optimal_result.sub_cluster_runtime, digits=4), " seconds")
    println("MFC Optimal Completion Edges Runtime:", round(mfc_optimal_result.completion_edges_runtime, digits=4), " seconds")
    println("MFC Optimal Total Runtime:           ", round(mfc_optimal_result.completion_runtime, digits=4), " seconds")

    println("\nMFC Simple Sub-cluster Runtime:      ", round(mfc_simple_result.sub_cluster_runtime, digits=4), " seconds")
    println("MFC Simple Completion Edges Runtime: ", round(mfc_simple_result.completion_edges_runtime, digits=4), " seconds")
    println("MFC Simple Total Runtime:            ", round(mfc_simple_result.completion_runtime, digits=4), " seconds")

    # Note that cookingndataset does not have ground truth labels
    # So can't really NMI, ARI scores

end

run_cooking_test("data/cooking.txt")