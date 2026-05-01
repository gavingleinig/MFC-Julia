using Plots
using DataFrames
using CSV



function load_matrix(gausin::Vector{Int64}, dim::Vector{Int64}, name::String)
    matrix_df = Vector{Vector{DataFrame}}()


    

    for dim_num in dim
        println("Testing dim: ",dim_num)
        dim_row = Vector{DataFrame}()
        for gausin_num in gausin
            println("Testing gausin: ", gausin_num)
            
            df = CSV.read(string(name,"_dim_",dim_num,"_gaus_",gausin_num,".csv"), DataFrame)
            push!(dim_row, df)
        end
        push!(matrix_df,dim_row)
    end

    return matrix_df
end




function plot_matrix_nmi(matrix_df::Vector{Vector{DataFrame}},gausin::Vector{Int64}, dim::Vector{Int64})
    
    final_plot = create_final_plot("NMI")
    plot_data(matrix_df, gausin, dim, final_plot,plot_dataframe_nmi)

    output_image = string("clustering_metrics_vs_sigma_dim_",length(dim),"_Gaus_",length(gausin),"_NMI",".png")
    savefig(final_plot, output_image)
    println("Plot saved to $output_image")

   



end


function plot_matrix_ari(matrix_df::Vector{Vector{DataFrame}},gausin::Vector{Int64}, dim::Vector{Int64})

    final_plot = create_final_plot("ARI")


    plot_data(matrix_df, gausin, dim, final_plot,plot_dataframe_ari)

    output_image = string("clustering_metrics_vs_sigma_dim_",length(dim),"_Gaus_",length(gausin),"_ARI",".png")
    savefig(final_plot, output_image)
    println("Plot saved to $output_image")

end

function plot_matrix_runtime(matrix_df::Vector{Vector{DataFrame}},gausin::Vector{Int64}, dim::Vector{Int64})

     l = @layout [grid(length(dim), length(gausin))
             b{0.05h}] 

    
    final_plot = plot(
        layout=l,
        xlabel ="",
        ylabel = "Runtime (sec)",
        legend = false,
        grid = true,
        xguidefontsize=4,
        yguidefontsize=4,
        xtickfontsize=4,
        ytickfontsize=4,
        titlefontsize=4

        )
    plot_data(matrix_df, gausin, dim, final_plot,plot_dataframe_runtime)

    output_image = string("clustering_metrics_vs_sigma_dim_",length(dim),"_Gaus_",length(gausin),"_ARI",".png")
    savefig(final_plot, output_image)
    println("Plot saved to $output_image")

end


function plot_matrix_gamma(matrix_df::Vector{Vector{DataFrame}},gausin::Vector{Int64}, dim::Vector{Int64})

     l = @layout [grid(length(dim), length(gausin))
             b{0.05h}] 

    
    final_plot = plot(
        layout=l,
        xlabel ="",
        ylabel = "ARI",
        legend = false,
        grid = true,
        xguidefontsize=4,
        yguidefontsize=4,
        xtickfontsize=4,
        ytickfontsize=4,
        titlefontsize=4,
        yticks=0:0.25:1,
        ylims=(0, 1),

        )
    plot_data(matrix_df, gausin, dim, final_plot,plot_dataframe_gamma)

    output_image = string("clustering_metrics_vs_sigma_dim_",length(dim),"_Gaus_",length(gausin),"_ARI",".png")
    savefig(final_plot, output_image)
    println("Plot saved to $output_image")

end

function matrix_table(matrix_df::Vector{Vector{DataFrame}},gausin::Vector{Int64}, dim::Vector{Int64})

     l = @layout [grid(length(dim), length(gausin))
             b{0.05h}] 

    
    final_plot = plot(
        layout=l,
        xlabel ="",
        ylabel = "ARI",
        legend = false,
        grid = true,
        xguidefontsize=4,
        yguidefontsize=4,
        xtickfontsize=4,
        ytickfontsize=4,
        titlefontsize=4,
        yticks=0:0.25:1,
        ylims=(0, 1),

        )
    plot_data(matrix_df, gausin, dim, final_plot,print_out_tabels)

    output_image = string("clustering_metrics_vs_sigma_dim_",length(dim),"_Gaus_",length(gausin),"_ARI",".png")
    savefig(final_plot, output_image)
    println("Plot saved to $output_image")

end



function plot_data(matrix_df::Vector{Vector{DataFrame}},gausin::Vector{Int64}, dim::Vector{Int64}, final_plot, plot_dataframe::F)where {F<:Function}


    for (i, dim_num) in enumerate(dim)
        for (j, gausin_num) in enumerate(gausin)
            label = ""
            gausin_label = ""

            if j == length(gausin)
                label = "Sigma"
            end

            if i == 1
                gausin_label = string("Gaussians: ",gausin_num,"\n")

            end
            if j == round(length(gausin)/2)
                gausin_label = string("Uniform Random Pontis d=",dim_num)
                if i == 1
                    gausin_label = string("Gaussians: ",gausin_num,"\nUniform Random Pontis d=",dim_num)
                end

                
                
            end
            
                
                    gausin_label = string("Gaussians: ",gausin_num,"\nd=",dim_num)
            

            println("Testing subplot: ", j+((i-1)*length(gausin)))
            plot_dataframe(matrix_df[i][j],final_plot,j+((i-1)*length(gausin)), label,gausin_label)
            
            
        end
        
        
        
    end

    # add ledger 
    plot!(final_plot, rand(0), subplot = length(gausin)*length(dim)+1, label = "K-Centering", marker = :circle,    legend = :inside,legendfontsize = 3, legend_column = 7, legendframestyle = :none, framestyle = :none, axis = false, grid = false)
    plot!(final_plot, rand(0), subplot = length(gausin)*length(dim)+1, label = "Optimal MST", marker = :square,    legend = :inside)
   plot!(final_plot, rand(0), subplot = length(gausin)*length(dim)+1, label = "Naive ST",    marker = :dtriangle, legend = :inside)
    plot!(final_plot, rand(0), subplot = length(gausin)*length(dim)+1, label = "MFC Approx",  marker = :diamond,   legend = :inside)
    plot!(final_plot, rand(0), subplot = length(gausin)*length(dim)+1, label = "MFC Optimal", marker = :star5,     legend = :inside)
    plot!(final_plot, rand(0), subplot = length(gausin)*length(dim)+1, label = "MFC Simple",  marker = :cross,     legend = :inside)
    #plot!(final_plot, rand(0), subplot = length(gausin)*length(dim)+1, label = "K-Means",     marker = :xcross,    legend = :inside)
    display(final_plot)


end

function create_final_plot(name::String)


    l = @layout [grid(length(dim), length(gausin))
             b{0.05h}] 

    
    final_plot = plot(
        layout=l,
        xlabel ="",
        ylabel = name,
        legend = false,
        grid = true,
        xguidefontsize=4,
        yguidefontsize=4,
        xtickfontsize=3,
        ytickfontsize=4,
        yticks=0:0.25:1,
        ylims=(0, 1),
        titlefontsize=4

        )

    return final_plot


end



function plot_dataframe_nmi(df::DataFrame,plots,j,label,gausin_label)
    sigmas = df.Sigma
    
    # NMI
    
    plot!(plots, sigmas,markersize = 2,  df.KC_NMI, marker=:circle, lw=1,  subplot=j, xlabel =label, title=gausin_label )
    plot!(plots, sigmas,markersize = 2,  df.MST_NMI, marker=:square, lw=1, subplot=j)
    plot!(plots, sigmas,markersize = 2,  df.Naive_NMI, marker=:dtriangle, lw=1, subplot=j)
    plot!(plots, sigmas,markersize = 2,  df.MFC_Approx_NMI, marker=:diamond, lw=1, subplot=j)
    plot!(plots, sigmas,markersize = 2,  df.MFC_Optimal_NMI, marker=:star5, lw=1, subplot=j)
    plot!(plots, sigmas,markersize = 2,  df.MFC_Simple_NMI, marker=:cross, lw=1, subplot=j)
    #plot!(plots, sigmas,markersize = 2,  df.KMeans_NMI, marker=:xcross, lw=1, subplot=j)


    # plot!(plots,p_ari, layout=(1, 2), size=(1000, 450), margin=5Plots.mm, subplot=j)
    


end

function plot_dataframe_ari(df::DataFrame,plots,j,label,gausin_label)
    sigmas = df.Sigma
    
    # ARI PLOT

    
    plot!(plots, sigmas,subplot=j,markersize = 2 , df.KC_ARI, marker=:circle, lw=.5, xlabel=label, title=gausin_label)
    plot!(plots, sigmas,subplot=j,markersize = 2 , df.MST_ARI, marker=:square, lw=.5)
    plot!(plots, sigmas,subplot=j,markersize = 2 , df.Naive_ARI, marker=:dtriangle, lw=.5)
    plot!(plots, sigmas,subplot=j,markersize = 2 , df.MFC_Approx_ARI, marker=:diamond, lw=.5)
    plot!(plots, sigmas,subplot=j,markersize = 2 , df.MFC_Optimal_ARI, marker=:star5, lw=.5)
    plot!(plots, sigmas,subplot=j,markersize = 2 , df.MFC_Simple_ARI, marker=:cross, lw=.5)
    #plot!(plots, sigmas,subplot=j,markersize = 2 , df.KMeans_ARI, marker=:xcross, lw=.5)



    # plot!(plots,p_ari, layout=(1, 2), size=(1000, 450), margin=5Plots.mm, subplot=j)
    


end


function plot_dataframe_gamma(df::DataFrame,plots,j,label,gausin_label)
    gamma = df.Gamma_Overlap
    
    # ARI PLOT

    
    # plot!(plots, gamma,subplot=j,markersize = 2 , df.KC_ARI, marker=:circle, lw=.5, xlabel=label, title=gausin_label)
    # plot!(plots, gamma,subplot=j,markersize = 2 , df.MST_ARI, marker=:square, lw=.5)
    #plot!(plots, sigmas,subplot=j,markersize = 2 , df.Naive_ARI, marker=:dtriangle, lw=.5)
    plot!(plots, gamma,subplot=j,markersize = 2 , df.MFC_Approx_NMI, marker=:diamond, lw=.5)
    # plot!(plots, gamma,subplot=j,markersize = 2 , df.MFC_Optimal_ARI, marker=:star5, lw=.5)
    #plot!(plots, sigmas,subplot=j,markersize = 2 , df.MFC_Simple_ARI, marker=:cross, lw=.5)
    #plot!(plots, sigmas,subplot=j,markersize = 2 , df.KMeans_ARI, marker=:xcross, lw=.5)



    # plot!(plots,p_ari, layout=(1, 2), size=(1000, 450), margin=5Plots.mm, subplot=j)
    


end


function plot_dataframe_runtime(df::DataFrame,plots,j,label,gausin_label)
    sigmas = df.Sigma
    
    # ARI PLOT

    
    plot!(plots, sigmas,subplot=j,markersize = 2 , df.KC_runtime, marker=:circle, lw=.5, xlabel=label, title=gausin_label)
    plot!(plots, sigmas,subplot=j,markersize = 2 , df.MST_runtime, marker=:square, lw=.5)
    plot!(plots, sigmas,subplot=j,markersize = 2 , df.Naive_runtime, marker=:dtriangle, lw=.5)
    plot!(plots, sigmas,subplot=j,markersize = 2 , df.MFC_Approx_runtime, marker=:diamond, lw=.5)
    plot!(plots, sigmas,subplot=j,markersize = 2 , df.MFC_Optimal_runtime, marker=:star5, lw=.5)
    plot!(plots, sigmas,subplot=j,markersize = 2 , df.MFC_Simple_runtime, marker=:cross, lw=.5)
    plot!(plots, sigmas,subplot=j,markersize = 2 , df.KMeans_runtime, marker=:xcross, lw=.5)



    # plot!(plots,p_ari, layout=(1, 2), size=(1000, 450), margin=5Plots.mm, subplot=j)
    


end


function print_out_tabels(df::DataFrame,plots,j,label,gausin_label)
    println(gausin_label)
    println(describe(df))

end

# ARI: plot_matrix_ari
# NMI: plot_matrix_nmi
# Tables: matrix_table
# Runtime: plot_matrix_runtime
#           NOTE: legacy Runtime is stored as runtime_[FILE NAME]
#                 New generated runtimes do not need to specify this
# Gamma vs. NMI: plot_matrix_gamma


gausin = [16,64]
# gausin = [16,32,64,128,256] # for full sweep

dim = [4,16,64]
# dim = [4,16,64,128] # for full sweep

output_file_name = "final_clustering_results" 
# output_file_name = "clustering_results" # for results used in paper
# output_file_name = "sigma_clustering_results" # for data with sigma to 2

# create graphs for cluster
df_matrix  = load_matrix(gausin,dim,string("data/Systhetic/",output_file_name))
plot_matrix_nmi(df_matrix,gausin,dim)



# LEGEACY Runtime example:
# gausin = [16,32,64,128,256]
# dim = [4,16,64,128]
# create graphs for runtime
# df_matrix  = load_matrix(gausin,dim,"runtime_clustering_results")
# matrix_table(df_matrix,gausin,dim)

