struct ClusteringResult{T<:Real}
    assignments::Vector{Int}
    runtime::T
end


"""
    closest_point(point, centroids::AbstractVector, dist_func::Function)

Find the index of the nearest centroid to a given 'point' using the provided distance function
"""
function closest_point(
    point,
    centroids::AbstractVector, 
    dist_func::Function
    )

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

"""
    k_centering(points::AbstractVector, k::Int, initial_index::Int, dist_func::Function)

Perform greedy 2-approx k-center clustering by iteratively picking the point furthest from all existing centroids. 
Returns a 'ClusteringResult' containing cluster assignments and runtime.
"""
function k_centering(
    points::AbstractVector, 
    k::Int,
    initial_index::Int,
    dist_func::F
) where {F<:Function}
    beginTime = time()

    if length(points) < k
        return nothing
    end
    
    if k <= 1
        return ClusteringResult(ones(Int, length(points)), 0.0)
    end

    
    centroids = Vector{eltype(points)}()
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
    return ClusteringResult(res, endTime - beginTime)
end



# function k_means()


# end

