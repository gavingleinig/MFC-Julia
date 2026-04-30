using MetricForestCompletion
using Random
using Clustering

function run_cooking_test(file_path::String, ground_truth_file_path::String)

    num_clusters = 20  # note that there are 20 ground truth clustering in the Cooking dataset
    dist_func = jaccard

    points = Vector{Int}[]
    ground_truth = Int[]
    edge_size_filter = 0 # Potentially filter out small sets
    
    open(file_path, "r") do f_data
        open(ground_truth_file_path, "r") do f_labels
            for (line, label_line) in zip(eachline(f_data), eachline(f_labels))
                parts = split(strip(line), ',')

                if isempty(parts) || parts == [""]
                    continue
                end
                
                s = sort!([parse(Int, p) for p in parts])
                
                if length(s) >= edge_size_filter
                    push!(points, s)
                    push!(ground_truth, parse(Int, strip(label_line)))
                end
            end
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

    initial_forest_complete = mfc_approx_result.cluster_edges
    gamma_overlap = compute_gamma_overlap(
        initial_forest_complete, 
        mst_complete, 
        clustering_result.assignments
    )

    naive_st_edges = naive_random_st(points, dist_func)
    naive_st_complete = convert_mst_to_weight(naive_st_edges)

        # Parameter sweeps for best thresholds
    t_sweep_mst = @elapsed best_result_mst, max_ari_mst = best_single_linkage_threshold_increment(mst_complete, num_points, ground_truth, 1.0)
    t_sweep_naive = @elapsed best_result_naive, max_ari_naive = best_single_linkage_threshold_increment(naive_st_complete, num_points, ground_truth, 1.0)

    t_sweep_approx = @elapsed best_result_approx, max_ari_approx = best_single_linkage_threshold_increment(all_ST_approx_edges, num_points, ground_truth, 1.0)
    t_sweep_optimal = @elapsed best_result_optimal, max_ari_optimal = best_single_linkage_threshold_increment(all_ST_optimal_edges, num_points, ground_truth, 1.0)
    t_sweep_simple = @elapsed best_result_simple, max_ari_simple = best_single_linkage_threshold_increment(all_ST_simple_edges, num_points, ground_truth, 1.0)


    # Cost ratio: How close are approx MFC trees to optimal ST?
    weight_mst = sum(e[3] for e in mst_complete)
    weight_approx = sum(e[3] for e in all_ST_approx_edges)
    weight_optimal = sum(e[3] for e in all_ST_optimal_edges)
    
    cost_ratio_approx = weight_approx / weight_mst
    cost_ratio_optimal = weight_optimal / weight_mst

    # What edges in  Optimal MST connect different ground truth cuisines?
    # How big / small are the distances are between cuisines?
    cross_cluster_edges = Tuple{Int, Int, Float64}[]
    for (u, v, w) in mst_complete
        if ground_truth[u] != ground_truth[v]
            push!(cross_cluster_edges, (u, v, w))
        end
    end
    
    sort!(cross_cluster_edges, by = x -> x[3])
    num_bridges = length(cross_cluster_edges)
    mean_bridge_weight = num_bridges > 0 ? sum(e[3] for e in cross_cluster_edges) / num_bridges : 0.0

    # What if we cut k longest edges? Is this any better?
    t_k_mst = @elapsed result_k_mst = single_linkage_k_clusters(mst_complete, num_points, num_clusters)
    t_k_approx = @elapsed result_k_approx = single_linkage_k_clusters(all_ST_approx_edges, num_points, num_clusters)
    t_k_optimal = @elapsed result_k_optimal = single_linkage_k_clusters(all_ST_optimal_edges, num_points, num_clusters)

    

    println("\n### Evaluation OLD ###")

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

    
    println("\n### Parameter Sweep Runtimes ###")
    println("MST Sweep Time:                      ", round(t_sweep_mst, digits=4), " seconds")
    println("Naive ST Sweep Time:                 ", round(t_sweep_naive, digits=4), " seconds")
    println("MFC Approx Sweep Time:               ", round(t_sweep_approx, digits=4), " seconds")
    println("MFC Optimal Sweep Time:              ", round(t_sweep_optimal, digits=4), " seconds")
    println("MFC Simple Sweep Time:               ", round(t_sweep_simple, digits=4), " seconds")

    println("\n### Evaluation ###")
    
    println("\nInitial K-Centering:")
    kc_ari = randindex(clustering_result.assignments, ground_truth)[1]
    kc_nmi = mutualinfo(clustering_result.assignments, ground_truth, normed=true)
    println("ARI: ", round(kc_ari, digits=4))
    println("NMI: ", round(kc_nmi, digits=4))

    println("y-overlap Bound:             ", round(gamma_overlap, digits=4))

    println("\nMST Single Linkage (Best Threshold):")
    mst_nmi = mutualinfo(best_result_mst.assignments, ground_truth, normed=true)
    println("ARI: ", round(max_ari_mst, digits=4))
    println("NMI: ", round(mst_nmi, digits=4))

    println("\nNaive ST Single Linkage (Best Threshold):")
    naive_nmi = mutualinfo(best_result_naive.assignments, ground_truth, normed=true)
    println("ARI: ", round(max_ari_naive, digits=4))
    println("NMI: ", round(naive_nmi, digits=4))

    println("\nMFC Approx Single Linkage (Best Threshold):")
    mfc_best_approx_nmi = mutualinfo(best_result_approx.assignments, ground_truth, normed=true)
    println("ARI: ", round(max_ari_approx, digits=4))
    println("NMI: ", round(mfc_best_approx_nmi, digits=4))

    println("\nMFC Optimal Single Linkage (Best Threshold):")
    mfc_best_optimal_nmi = mutualinfo(best_result_optimal.assignments, ground_truth, normed=true)
    println("ARI: ", round(max_ari_optimal, digits=4))
    println("NMI: ", round(mfc_best_optimal_nmi, digits=4))

    println("\nMFC Simple Single Linkage (Best Threshold):")
    mfc_best_simple_nmi = mutualinfo(best_result_simple.assignments, ground_truth, normed=true)
    println("ARI: ", round(max_ari_simple, digits=4))
    println("NMI: ", round(mfc_best_simple_nmi, digits=4))

    println("\n### Evaluation Pt 2###")
    
    println("\nCost Rations")
    println("MFC Approx Cost Ratio:       ", round(cost_ratio_approx, digits=4))
    println("MFC Optimal Cost Ratio:      ", round(cost_ratio_optimal, digits=4))
    
    println("\nMST Cluster Analysis:")
    println("Total cross-cuisine edges:   ", num_bridges)
    if num_bridges > 0
        println("Shortest cross-cuisine edge: ", round(cross_cluster_edges[1][3], digits=4))
        println("Average cross-cuisine edge:  ", round(mean_bridge_weight, digits=4))
    end

    println("\nSingle Linkage w/ cut 19 Longest Edges")
    println("Optimal MST ARI: ", round(randindex(result_k_mst.assignments, ground_truth)[1], digits=4), "NMI: ", round(mutualinfo(result_k_mst.assignments, ground_truth, normed=true), digits=4))
    println("MFC Approx ARI: ", round(randindex(result_k_approx.assignments, ground_truth)[1], digits=4), "NMI: ", round(mutualinfo(result_k_approx.assignments, ground_truth, normed=true), digits=4))


end

run_cooking_test("data/Cooking/cooking.txt","data/Cooking/cooking_labels.txt")