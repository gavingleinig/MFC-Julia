# k_means
# k_centering
# perhaps use Julia's Clustering.jl?

struct ClusteringResult{T<:AbstractFloat}
    assignments::Vector{Int}
    runtime::T
end

# V <: AbstractVector: Represents the type of a single point (e.g., Vector{Float64})
# F <: Function: Specializes the compiler for the exact distance function provided
function k_centering(
    points::AbstractVector{V}, 
    k::Int, 
    dist_func::F
) where {V<:AbstractVector, F<:Function}
    
    # TODO
    
    return ClusteringResult{Float64}(Int[], 0.0) 
end