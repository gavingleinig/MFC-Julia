using MetricForestCompletion
using Random
using Clustering
using Statistics
using DataFrames
using CSV

# Runs every procedure `num_trials` times on freshly generated Gaussian data
# and reports the mean (and std) of ARI, NMI and runtime for each method.
function run_gaussian_trials(;
    num_trials::Int = 30,
    base_seed::Int = 1000,
    precision::Type{<:AbstractFloat} = Float64,
    dim::Int = 8,
    num_gauss::Int = 16,
    num_points::Int = 20000,
    num_clusters::Int = num_gauss,          # defaults to one cluster per Gaussian
    mean_range::Tuple{<:Real, <:Real} = (-5, 5),
    sigma_range::Tuple{<:Real, <:Real} = (0.5, 0.8),
    threshold_step::Real = 0.01,            # increment for the threshold sweep
    dist_func = euclidean,
    output_dir::String = "data/Systhetic/",
    csv_label::String = "gaussian_trials"
)
    PrecisionType = precision
    mean_range = (PrecisionType(mean_range[1]), PrecisionType(mean_range[2]))
    sigma_range = (PrecisionType(sigma_range[1]), PrecisionType(sigma_range[2]))

    methods = ["KMeans", "KC", "MST", "Naive", "MFC_Approx", "MFC_Optimal", "MFC_Simple", "MFC_Multi"]
    metrics = ["ARI", "NMI", "runtime"]

    cols = Dict{String, Vector{Float64}}()
    for m in methods, metric in metrics
        cols["$(m)_$(metric)"] = Float64[]
    end
    gamma_overlaps = Float64[]
    seeds = Int[]

    function record!(m, assignments, ari, runtime, ground_truth)
        push!(cols["$(m)_ARI"], ari)
        push!(cols["$(m)_NMI"], mutualinfo(assignments, ground_truth, normed=true))
        push!(cols["$(m)_runtime"], runtime)
    end

    for trial in 1:num_trials
        seed = base_seed + trial
        Random.seed!(seed)  # seeds any global RNG use inside the library (e.g. kmeans)
        println("Trial $trial / $num_trials (seed = $seed)")
        trial_start = time()
        push!(seeds, seed)

        # A fresh random dataset each trial (new means/sigmas/points)
        points, ground_truth = generate_gaussians(
            PrecisionType;
            dim = dim,
            num_centers = num_gauss,
            points_per_center = num_points ÷ num_gauss,
            mean_range = mean_range,
            sigma_range = sigma_range,
            rng = Xoshiro(seed)
        )
        n_pts = length(points)

        # Initial clustering and full k-means
        k_centering_runtime = @elapsed clustering_result = k_centering(points, num_clusters, 1, dist_func)
        points_matrix = stack(points)
        kmeans_runtime = @elapsed kmeans_result = kmeans(points_matrix, num_clusters)

        # MFC variants
        mfc_approx_runtime  = @elapsed mfc_approx_result  = metric_forest_completion_approx(points, num_clusters, clustering_result.assignments, dist_func)
        mfc_optimal_runtime = @elapsed mfc_optimal_result = metric_forest_completion_optimal(points, num_clusters, clustering_result.assignments, dist_func)
        mfc_simple_runtime  = @elapsed mfc_simple_result  = metric_forest_completion_simple(points, num_clusters, clustering_result.assignments, dist_func)
        mfc_multi_runtime   = @elapsed mfc_multi_result   = mulit_reps_metric_forest_completion(points, num_clusters, clustering_result.assignments, dist_func)

        approx_edges  = vcat(mfc_approx_result.cluster_edges,  mfc_approx_result.completion_edges)
        optimal_edges = vcat(mfc_optimal_result.cluster_edges, mfc_optimal_result.completion_edges)
        simple_edges  = vcat(mfc_simple_result.cluster_edges,  mfc_simple_result.completion_edges)
        multi_edges   = vcat(mfc_multi_result.cluster_edges,   mfc_multi_result.completion_edges)

        # Baselines: optimal MST and naive random ST
        mst_runtime   = @elapsed mst_complete   = convert_mst_to_weight(mst_implicit(points, dist_func))
        naive_runtime = @elapsed naive_complete = convert_mst_to_weight(naive_random_st(points, dist_func))

        push!(gamma_overlaps, Float64(compute_gamma_overlap(
            mfc_approx_result.cluster_edges, mst_complete, clustering_result.assignments)))

        # Best-threshold sweeps (best ARI vs. ground truth)
        best_mst,     mst_ari     = best_single_linkage_threshold_increment(mst_complete,   n_pts, ground_truth, threshold_step)
        best_naive,   naive_ari   = best_single_linkage_threshold_increment(naive_complete, n_pts, ground_truth, threshold_step)
        best_approx,  approx_ari  = best_single_linkage_threshold_increment(approx_edges,   n_pts, ground_truth, threshold_step)
        best_optimal, optimal_ari = best_single_linkage_threshold_increment(optimal_edges,  n_pts, ground_truth, threshold_step)
        best_simple,  simple_ari  = best_single_linkage_threshold_increment(simple_edges,   n_pts, ground_truth, threshold_step)
        best_multi,   multi_ari   = best_single_linkage_threshold_increment(multi_edges,    n_pts, ground_truth, threshold_step)

        # Record (MFC runtimes include the initial k-centering, as in evaluation.jl)
        record!("KMeans", kmeans_result.assignments, randindex(kmeans_result.assignments, ground_truth)[1], kmeans_runtime, ground_truth)
        record!("KC", clustering_result.assignments, randindex(clustering_result.assignments, ground_truth)[1], k_centering_runtime, ground_truth)
        record!("MST", best_mst.assignments, mst_ari, mst_runtime, ground_truth)
        record!("Naive", best_naive.assignments, naive_ari, naive_runtime, ground_truth)
        record!("MFC_Approx", best_approx.assignments, approx_ari, mfc_approx_runtime + k_centering_runtime, ground_truth)
        record!("MFC_Optimal", best_optimal.assignments, optimal_ari, mfc_optimal_runtime + k_centering_runtime, ground_truth)
        record!("MFC_Simple", best_simple.assignments, simple_ari, mfc_simple_runtime + k_centering_runtime, ground_truth)
        record!("MFC_Multi", best_multi.assignments, multi_ari, mfc_multi_runtime + k_centering_runtime, ground_truth)

        println("   Completed in $(round(time() - trial_start, digits=2)) sec")
    end

    # Per-trial table
    trial_df = DataFrame(
        [:Trial => collect(1:num_trials); :Seed => seeds; :Gamma_Overlap => gamma_overlaps;
         [Symbol(k) => cols[k] for m in methods for k in ["$(m)_ARI", "$(m)_NMI", "$(m)_runtime"]]]
    )

    # Summary table: mean (and std) over trials
    summary_df = DataFrame(
        Method       = methods,
        Mean_ARI     = [mean(cols["$(m)_ARI"]) for m in methods],
        Std_ARI      = [std(cols["$(m)_ARI"]) for m in methods],
        Mean_NMI     = [mean(cols["$(m)_NMI"]) for m in methods],
        Std_NMI      = [std(cols["$(m)_NMI"]) for m in methods],
        Mean_runtime = [mean(cols["$(m)_runtime"]) for m in methods],
    )

    println("\n### Gaussian Trials: mean over $num_trials runs ###")
    println("Mean gamma overlap: ", round(mean(gamma_overlaps), digits=4))
    println(summary_df)

    mkpath(output_dir)
    CSV.write(joinpath(output_dir, csv_label * "_per_trial.csv"), trial_df)
    CSV.write(joinpath(output_dir, csv_label * "_summary.csv"), summary_df)
    println("\nResults saved to $(output_dir)")

    return trial_df, summary_df
end

trial_df, summary_df = run_gaussian_trials()
