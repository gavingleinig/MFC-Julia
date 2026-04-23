using Plots
using DataFrames
using CSV

function ploting_data(file_path::String)
    
    df = CSV.read(file_path, DataFrame)
    println(df.Sigma)
    p_ARI = plot(df.Sigma, [df.KC_ARI,df.MST_ARI,df.Naive_ARI,df.MFC_Approx_ARI,df.MFC_Optimal_ARI,df.MFC_Simple_ARI], title = "Adjusted Rand Index (ARI) vs Sigma", label=["K-Centering" "Optimal MST" "Naive ST" "MFC Approx" "MFC Optimal" "MFC Simple"], titlefont = font(8))
    p_NMI = plot(df.Sigma, [df.KC_NMI,df.MST_NMI,df.Naive_NMI,df.MFC_Approx_NMI,df.MFC_Optimal_NMI,df.MFC_Simple_NMI], label=["K-Centering" "Optimal MST" "Naive ST" "MFC Approx" "MFC Optimal" "MFC Simple"], title = "Normalized Mutual Info (NMI) vs Sigma",titlefont = font(8))
    plot(p_ARI, p_NMI)


end

ploting_data("clustering_results.csv")