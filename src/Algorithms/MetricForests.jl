# src/Algorithms/MetricForests.jl

struct MFCResult{T<:AbstractFloat}
    cluster_edges::Vector{Tuple{Int, Int}}      
    completion_edges::Vector{Tuple{Int, Int}}   
    sub_cluster_runtime::T
    completion_edges_runtime::T
    completion_runtime::T
end

function metric_forest_completion(
    points::AbstractVector{V}, 
    cluster_count::Int, 
    cluster_assignments::Vector{Int}, 
    dist_func::F
) where {V<:AbstractVector, F<:Function}
    
    # TODO
    
    return MFCResult{Float64}(Tuple{Int, Int}[], Tuple{Int, Int}[], 0.0, 0.0, 0.0) 
end