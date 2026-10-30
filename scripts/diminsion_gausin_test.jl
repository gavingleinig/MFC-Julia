using MetricForestCompletion
using Random
using Clustering
using DataFrames
using CSV
 
include("evaluation.jl")    # defines run_gaussian_sigma_sweep
include("rings_trials.jl")  # defines run_rings_trials
 
 
function dimension_gaussian_test()
    # true -> run the rings trials, false -> run the Gaussian sigma sweep
    rings = true
 
    # change number of gaussians (or rings)
    gaussians = [3, 4, 5]
    # gaussians = [16, 32, 64, 128, 256]  # for full sweep
 
    # change number of dimensions
    # dims = [2, 3, 4]
    dims = [4, 16, 64, 128]  # for full sweep
 
    # change number of points
    points = 2000
    # points = 20000  # for full test
 
    # change mean range (Gaussian sweep only)
    mean_range = (-5, 5)
 
    # change name
    output_name = "final_clustering_results"
 
    # change sigma range (Gaussian sweep only)
    sigmas_to_test = 0.1:0.2:1.5
 
    for dim_num in dims
        println("Testing dim: ", dim_num)
 
        for gaussian_num in gaussians
            println("Testing gaussians: ", gaussian_num)
 
            # specify dim and gaussians
            output_file_name_full = string(output_name, "_dim_", dim_num, "_gaus_", gaussian_num, ".csv")
 
            if rings
                run_rings_trials(
                    dim = dim_num,
                    num_rings = gaussian_num,
                    points_per_ring = points ÷ gaussian_num,
                    csv_label = replace(output_file_name_full, ".csv" => "")
                )
            else
                run_gaussian_sigma_sweep(
                    sigmas_to_test, dim_num, gaussian_num,
                    output_file_name_full, points, mean_range
                )
            end
        end
    end
end
 
 
dimension_gaussian_test()