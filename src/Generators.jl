# Synthetic Data Generation
# Replace gaussian_point_gen.cpp and unifrom.cpp
using Random
using Distributions

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
    
    total_points = num_centers * points_per_center
    points = Vector{Vector{T}}(undef, total_points)
    labels = Vector{Int}(undef, total_points)
    # gerneate the random number
    g = Xoshiro(123)


    
    idx = 1
    for cluster_id in 1:num_centers
        # TODO: Setup standard normal distributions for this cluster_id w/ mean_range and sigma_range
        
        # First randimoly generate Mean for cluster
        cluster_mean = (rand(g) * (mean_range[2] - mean_range[1] + 1)) + mean_range[1]
        # Next randimly genearte sigma for cluster
        cluster_sigma = (rand(g) * (sigma_range[2] - sigma_range[1] + 1)) + sigma_range[1]
        # Finally create nomral distibution
        d = Normal(cluster_mean, cluster_sigma)

        for p in 1:points_per_center
            # TODO: Generate individual point 
            points[idx] = rand(d,dim)
            
            labels[idx] = cluster_id 
            idx += 1
        end
    end
    
    return points, labels
end