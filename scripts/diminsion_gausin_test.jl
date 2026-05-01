using MetricForestCompletion
using Random
using Clustering
using DataFrames
using CSV
include("evaluation.jl")



function diminsion_gausian_test()
    gausin = [16,64,256]
    dim = [4,16,32]

    for dim_num in dim
        println("Testing dim: ",dim_num)
        
        for gausin_num in gausin
            println("Testing gausin: ", gausin_num)
            sigmas_to_test = 0.1:0.2:2
            df_results = run_gaussian_sigma_sweep(sigmas_to_test,dim_num,gausin_num,string("sigma_clustering_results_dim_",dim_num,"_gaus_",gausin_num,".csv"))
        end
    end




end


diminsion_gausian_test()