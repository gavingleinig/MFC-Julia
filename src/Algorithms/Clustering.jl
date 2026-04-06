# k_means
# k_centering
# perhaps use Julia's Clustering.jl? Not a bad idea I will look into it
# we can further use diffrent clustering algorithems to see if we increase effincancy
using()

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
    beginTime = time()

    if size(points) < k
        return
    end
    
    if k <= 1
        return ClusteringResult{Float64}(points, 0) 
    end

    
    





    # TODO



    endTime = time()
    return ClusteringResult{Float64}(Int[], beginTime - endTime) 
end



function k_means()


end