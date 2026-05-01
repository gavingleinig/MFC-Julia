using MetricForestCompletion
using Random
using Clustering
using DataFrames
using CSV
include("evaluation.jl")



function diminsion_gausian_test()
    # change number of gausians

    
    gausin = [16,64]
    # gausin = [16,32,64,128,256] # for full sweep

    # change numer of diminsions
    dim = [4,16,64]
    # dim = [4,16,64,128] # for full sweep

    # change number of pionts:
    points = 2000 
    # points = 20000 #for full test
    
    # change mean range
    mean_range = (-5,5)
    # change name
    output_name = string("final_clustering_results")
    # Change sigma range
    sigmas_to_test = 0.1:0.2:1.5


    for dim_num in dim
        println("Testing dim: ",dim_num)
        
        for gausin_num in gausin
            println("Testing gausin: ", gausin_num)
            # specify dim and gausin
            output_file_name_full = string(output_name,"_dim_",dim_num,"_gaus_",gausin_num,".csv")


            df_results = run_gaussian_sigma_sweep(sigmas_to_test,dim_num,gausin_num,output_file_name_full,points,mean_range)
        end
    end




end


diminsion_gausian_test()