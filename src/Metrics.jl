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

function jaccard(a::AbstractVector{Int}, b::AbstractVector{Int})
    union_set = union(a, b)
    union_size = length(union_set)
    
    if union_size == 0
        return 1.0
    end
    
    intersection_size = length(intersect(a, b))
    return 1.0 - (intersection_size / union_size)
end

function hamming(a::AbstractString, b::AbstractString)
    if length(a) != length(b)
        throw(ArgumentError("Sequences not same length"))
    end
    
    distance = 0.0
    for (char_a, char_b) in zip(a, b)
        if char_a != char_b
            distance += 1.0
        end
    end
    
    return distance
end