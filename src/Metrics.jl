# Contains distance logic, 
# Include custom distances like edit_distance.cpp
# Perhaps wrap other metric from 'Distances.jl'

# Calculates the Euclidean distance between two points
function euclidean(a::AbstractVector{T}, b::AbstractVector{T}) where {T<:Real}
    # TODO
    i = 1
    distance = 0

    while  i <= length(a)
        distance += (a[i] - b[i])^2
        i+=1
    end

    distance = sqrt(distance)
    # return zero(T) 
    return distance
end

function jaccard(a::Set{Int}, b::Set{Int})
    union_set = union(a, b)
    union_size = length(union_set)
    
    if union_size == 0
        return 1.0
    end
    
    intersection_size = length(intersect(a, b))
    return 1.0 - (intersection_size / union_size)
end