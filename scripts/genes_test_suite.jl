using MetricForestCompletion
using Random
using Clustering
using Statistics

function run_genes_test_suite(file_path::String, ground_truth_file_path::String; num_runs::Int=5)
    dist_func = hamming

    taxonomies = [
        (2, "Phylum", 3),
        (3, "Class", 20),
        (4, "Order", 38),
        (5, "Family", 67),
        (6, "Genus", 130),
        (7, "Species", 148)
    ]

    base_seeds = [1862, 1876, 1939, 2025, 2027]
    seeds = base_seeds[1:min(num_runs, length(base_seeds))] 

    for (target_level, level_name, num_clusters) in taxonomies

        points = String[]
        ground_truth = Int[]
        edge_size_filter = 0 # Potentially filter out small sets
        
        open(file_path, "r") do f_data
            open(ground_truth_file_path, "r") do f_labels
                # Skip CSV header line
                readline(f_labels) 
                
                for (line, label_line) in zip(eachline(f_data), eachline(f_labels))
                    push!(points, strip(line))
                    
                    label_parts = split(strip(label_line), ',')
                    push!(ground_truth, parse(Int, label_parts[target_level]))
                end
            end
        end
        num_points = length(points)
        println("Loaded $(num_points) sets from genes dataset for $(level_name).")

        println("Computing Exact MST...")
        t_mst_compute = @elapsed mst_edges = mst_implicit(points, dist_func)
        mst_complete = convert_mst_to_weight(mst_edges)

        # Accuracy metric arrays
        ari_kc, nmi_kc = Float64[], Float64[]
        ari_mst, nmi_mst = Float64[], Float64[]
        ari_naive, nmi_naive = Float64[], Float64[]
        ari_mfc_approx, nmi_mfc_approx = Float64[], Float64[]
        ari_mfc_optimal, nmi_mfc_optimal = Float64[], Float64[]
        ari_mfc_simple, nmi_mfc_simple = Float64[], Float64[]

        # Runtime arrays for averaging
        rt_kc = Float64[]
        
        rt_mfc_approx_sub = Float64[]; rt_mfc_approx_comp = Float64[]; rt_mfc_approx_tot = Float64[]
        rt_mfc_optimal_sub = Float64[]; rt_mfc_optimal_comp = Float64[]; rt_mfc_optimal_tot = Float64[]
        rt_mfc_simple_sub = Float64[]; rt_mfc_simple_comp = Float64[]; rt_mfc_simple_tot = Float64[]
        
        rt_sweep_mst = Float64[]
        rt_sweep_naive = Float64[]
        rt_sweep_approx = Float64[]
        rt_sweep_optimal = Float64[]
        rt_sweep_simple = Float64[]

        for (run_idx, seed_val) in enumerate(seeds)
            # Set the random seed to vary start points for k_centering and naive ST
            Random.seed!(seed_val)
            println("Executing Run $(run_idx)/$(num_runs) (Seed: $(seed_val))")

            clustering_result = k_centering(points, num_clusters, rand(1:num_points), dist_func)
            push!(rt_kc, clustering_result.runtime)

            mfc_approx_result = metric_forest_completion_approx(points, num_clusters, clustering_result.assignments, dist_func)
            push!(rt_mfc_approx_sub, mfc_approx_result.sub_cluster_runtime)
            push!(rt_mfc_approx_comp, mfc_approx_result.completion_edges_runtime)
            push!(rt_mfc_approx_tot, mfc_approx_result.completion_runtime)

            mfc_optimal_result = metric_forest_completion_optimal(points, num_clusters, clustering_result.assignments, dist_func)
            push!(rt_mfc_optimal_sub, mfc_optimal_result.sub_cluster_runtime)
            push!(rt_mfc_optimal_comp, mfc_optimal_result.completion_edges_runtime)
            push!(rt_mfc_optimal_tot, mfc_optimal_result.completion_runtime)

            mfc_simple_result = metric_forest_completion_simple(points, num_clusters, clustering_result.assignments, dist_func)
            push!(rt_mfc_simple_sub, mfc_simple_result.sub_cluster_runtime)
            push!(rt_mfc_simple_comp, mfc_simple_result.completion_edges_runtime)
            push!(rt_mfc_simple_tot, mfc_simple_result.completion_runtime)

            all_ST_approx_edges  = vcat(mfc_approx_result.cluster_edges, mfc_approx_result.completion_edges)
            all_ST_optimal_edges = vcat(mfc_optimal_result.cluster_edges, mfc_optimal_result.completion_edges)
            all_ST_simple_edges  = vcat(mfc_simple_result.cluster_edges, mfc_simple_result.completion_edges)

            naive_st_edges = naive_random_st(points, dist_func)
            naive_st_complete = convert_mst_to_weight(naive_st_edges)

            # Parameter sweeps with elapsed times
            t_sweep_mst = @elapsed best_result_mst, max_ari_mst = best_single_linkage_threshold_increment(mst_complete, num_points, ground_truth, 1.0)
            push!(rt_sweep_mst, t_sweep_mst)

            t_sweep_naive = @elapsed best_result_naive, max_ari_naive = best_single_linkage_threshold_increment(naive_st_complete, num_points, ground_truth, 1.0)
            push!(rt_sweep_naive, t_sweep_naive)

            t_sweep_approx = @elapsed best_result_approx, max_ari_approx = best_single_linkage_threshold_increment(all_ST_approx_edges, num_points, ground_truth, 1.0)
            push!(rt_sweep_approx, t_sweep_approx)

            t_sweep_optimal = @elapsed best_result_optimal, max_ari_optimal = best_single_linkage_threshold_increment(all_ST_optimal_edges, num_points, ground_truth, 1.0)
            push!(rt_sweep_optimal, t_sweep_optimal)

            t_sweep_simple = @elapsed best_result_simple, max_ari_simple = best_single_linkage_threshold_increment(all_ST_simple_edges, num_points, ground_truth, 1.0)
            push!(rt_sweep_simple, t_sweep_simple)

            # Record metrics
            push!(ari_kc, randindex(clustering_result.assignments, ground_truth)[1])
            push!(nmi_kc, mutualinfo(clustering_result.assignments, ground_truth, normed=true))

            push!(ari_mst, max_ari_mst)
            push!(nmi_mst, mutualinfo(best_result_mst.assignments, ground_truth, normed=true))

            push!(ari_naive, max_ari_naive)
            push!(nmi_naive, mutualinfo(best_result_naive.assignments, ground_truth, normed=true))

            push!(ari_mfc_approx, max_ari_approx)
            push!(nmi_mfc_approx, mutualinfo(best_result_approx.assignments, ground_truth, normed=true))

            push!(ari_mfc_optimal, max_ari_optimal)
            push!(nmi_mfc_optimal, mutualinfo(best_result_optimal.assignments, ground_truth, normed=true))

            push!(ari_mfc_simple, max_ari_simple)
            push!(nmi_mfc_simple, mutualinfo(best_result_simple.assignments, ground_truth, normed=true))
        end

        println("\n==============================================")
        println(" Results over $(num_runs) runs (Mean) for $(level_name)")
        println("==============================================\n")
        
        println("### Clustering Metrics ###")
        println("\nInitial K-Centering:")
        println("  ARI: ", round(mean(ari_kc), digits=4))
        println("  NMI: ", round(mean(nmi_kc), digits=4))

        println("\nMST Single Linkage (Best Threshold):")
        println("  ARI: ", round(mean(ari_mst), digits=4))
        println("  NMI: ", round(mean(nmi_mst), digits=4))

        println("\nNaive ST Single Linkage (Best Threshold):")
        println("  ARI: ", round(mean(ari_naive), digits=4))
        println("  NMI: ", round(mean(nmi_naive), digits=4))

        println("\nMFC Approx Single Linkage (Best Threshold):")
        println("  ARI: ", round(mean(ari_mfc_approx), digits=4))
        println("  NMI: ", round(mean(nmi_mfc_approx), digits=4))

        println("\nMFC Optimal Single Linkage (Best Threshold):")
        println("  ARI: ", round(mean(ari_mfc_optimal), digits=4))
        println("  NMI: ", round(mean(nmi_mfc_optimal), digits=4))

        println("\nMFC Simple Single Linkage (Best Threshold):")
        println("  ARI: ", round(mean(ari_mfc_simple), digits=4))
        println("  NMI: ", round(mean(nmi_mfc_simple), digits=4))

        println("\n### Average Runtimes (Seconds) ###\n")
        println("Optimal MST Computation (Computed Once): ", round(t_mst_compute, digits=4))

        println("\n--- MFC Runtimes ---")
        println("Initial Clustering:                  ", round(mean(rt_kc), digits=4))
        
        println("\nMFC Approx Sub-cluster Runtime:      ", round(mean(rt_mfc_approx_sub), digits=4))
        println("MFC Approx Completion Edges Runtime: ", round(mean(rt_mfc_approx_comp), digits=4))
        println("MFC Approx Total Runtime:            ", round(mean(rt_mfc_approx_tot), digits=4))

        println("\nMFC Optimal Sub-cluster Runtime:     ", round(mean(rt_mfc_optimal_sub), digits=4))
        println("MFC Optimal Completion Edges Runtime:", round(mean(rt_mfc_optimal_comp), digits=4))
        println("MFC Optimal Total Runtime:           ", round(mean(rt_mfc_optimal_tot), digits=4))

        println("\nMFC Simple Sub-cluster Runtime:      ", round(mean(rt_mfc_simple_sub), digits=4))
        println("MFC Simple Completion Edges Runtime: ", round(mean(rt_mfc_simple_comp), digits=4))
        println("MFC Simple Total Runtime:            ", round(mean(rt_mfc_simple_tot), digits=4))

        println("\n--- Parameter Sweep Runtimes ---")
        println("MST Sweep Time:                      ", round(mean(rt_sweep_mst), digits=4))
        println("Naive ST Sweep Time:                 ", round(mean(rt_sweep_naive), digits=4))
        println("MFC Approx Sweep Time:               ", round(mean(rt_sweep_approx), digits=4))
        println("MFC Optimal Sweep Time:              ", round(mean(rt_sweep_optimal), digits=4))
        println("MFC Simple Sweep Time:               ", round(mean(rt_sweep_simple), digits=4))
    end
    
    println("\nTest Suite Completed.")
end

run_genes_test_suite("data/Genes/gg_13_5_ssualign_filtered_3000.txt", "data/Genes/gg_13_5_labels_all_3000.csv", num_runs=5)