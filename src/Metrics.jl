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