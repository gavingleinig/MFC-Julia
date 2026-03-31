# Synthetic Data Generation
# Replace gaussian_point_gen.cpp and unifrom.cpp

# Type{T} argument allows users to call generate_gaussians(Float32, dim=...) to get 32-bit floats
# Unsure if this is doing too much.
function generate_gaussians(
    ::Type{T}; 
    dim::Int, 
    num_centers::Int, 
    points_per_center::Int, 
    mean_range::Tuple{T, T}, 
    sigma_range::Tuple{T, T}
) where {T<:AbstractFloat}
    
    # TODO
    
    return Vector{Vector{T}}() 
end