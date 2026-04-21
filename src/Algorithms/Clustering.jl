# k_means
# k_centering
# perhaps use Julia's Clustering.jl? Not a bad idea I will look into it
# we can further use diffrent clustering algorithems to see if we increase effincancy
#using()

struct ClusteringResult{T<:AbstractFloat}
    assignments::Vector{Int}
    runtime::T
end


# helper functions

function closest_point(
    point::AbstractVector,
    centroids::AbstractVector{V}, 
    dist_func::Function
    )where {V<:AbstractVector}
    #find the closet point from the given point to a list of points
    res = 1
    dist = dist_func(point, centroids[1])
    
    for i in 2:length(centroids)
        cur_dist = dist_func(point, centroids[i])

        if cur_dist < dist 
            dist = cur_dist
            res = i
        end
    end
    return res
end

# V <: AbstractVector: Represents the type of a single point (e.g., Vector{Float64})
# F <: Function: Specializes the compiler for the exact distance function provided
# furtheset point function is never used so I am not including it
function k_centering(
    points::AbstractVector{V}, 
    k::Int,
    initial_index::Int,
    dist_func::F
) where {V<:AbstractVector, F<:Function}
    beginTime = time()

    if length(points) < k
        return nothing
    end
    
    if k <= 1
        return ClusteringResult{Float64}(ones(Int, length(points)), 0.0)
    end

    
    centroids = Vector{V}()
    sizehint!(centroids, k)
    push!(centroids, points[initial_index])

    cur_distances = Float64[]
    sizehint!(cur_distances, length(points))

    # calc initial distances of all points to first centroid
    for i in 1:length(points)
        push!(cur_distances, dist_func(points[i], centroids[1]))
    end

    # find other centroids
    while length(centroids) < k
        # find point with max distance to its nearest centroid; this is new centroid
        _, new_index = findmax(cur_distances)
        push!(centroids, points[new_index])

        # Update min distances to the new centroid
        for i in 1:length(points)
            dist = dist_func(points[i], centroids[end])
            if dist < cur_distances[i]
                cur_distances[i] = dist
            end
        end
    end

    res = Int[]
    sizehint!(res, length(points))

    for i in 1:length(points)
        push!(res, closest_point(points[i], centroids, dist_func))
    end

    endTime = time()
    return ClusteringResult{Float64}(res, endTime - beginTime)
end



# function k_means()


# end

