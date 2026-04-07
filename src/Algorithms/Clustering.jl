# k_means
# k_centering
# perhaps use Julia's Clustering.jl? Not a bad idea I will look into it
# we can further use diffrent clustering algorithems to see if we increase effincancy
using()

struct ClusteringResult{T<:AbstractFloat}
    assignments::Vector{Int}
    runtime::T
end


# helper functions

function closest_point(
    point, 
    points::V, 
    dist_func::F
    )where {V<:AbstractVector, F<:Function} 
    #find the closet point from the given point to a list of points
    res = 0
    dist = dist_func(point, points[0])
    i = 1 

    while i < length(V)
        cur_dist = dist_func(point, points[i])

        if cur_dist > dist 
            dist = cur_dist
            res = i
        end

        i++
    end

    return res
end

# V <: AbstractVector: Represents the type of a single point (e.g., Vector{Float64})
# F <: Function: Specializes the compiler for the exact distance function provided
# furtheset point function is never used so I am not including it
function k_centering(
    points::AbstractVector{V}, 
    k::Int,
    inital_index::Int,
    dist_func::F
) where {V<:AbstractVector, F<:Function}
    beginTime = time()

    if length(points) < k
        return nothing
    end
    
    if k <= 1
        return ClusteringResult{Float64}([length(points), 0], 0) 
    end

    
    centroids = AbstractVector{V}()
    sizehint!(centroids, num_clusters)
    push!(centroids, points[inital_index])


    cur_distances = AbstractVector()
    sizehint!(cur_distances, length(points))


    second_index = 0

    max_dist = dist_func(points[0], centroids[0])

    i = 0

    while i < length(points)
        dist = dist_func(points[i], centroids[0])
        if dist > max_dist
            second_index = i
            max_dist = dist
        end

        push!(cur_distances, dist)

        i++
    end

    push!(centroids, points[second_idex])


    while length(centroids) < k
        new_index = 0
        max_dist = minimum!([dist_func(points[0], centroids[length(centroids) - 1]), cur_distance[0]])
        cur_distance[0] = max_dist
        i  = 1
        while i < length(points)
            dist = dist_func(points[0], centroids[length(centroids) - 1])
            if dist < cur_distance[i]
                cur_distances[i] = dist
            end

            if cur_distances[i] > max_dist
                max_dist = cur_distances[i]
                new_index = i
            end

            i++
        end

        push!(centroids, points[new_index])
    end

    assignments = Vector{Int}()
    sizehint!(res, length(points))

    i = 0

    while i < length(points)
        push!(res, closet_point(points[i], centroids))
    end


    endTime = time()
    return ClusteringResult{Float64}(Int[], beginTime - endTime) 
end



function k_means()


end