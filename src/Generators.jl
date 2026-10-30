# Synthetic Data Generation
# Replace gaussian_point_gen.cpp and unifrom.cpp
using Random
using Distributions
using LinearAlgebra

# Type{T} argument allows users to call generate_gaussians(Float32, dim=...) to get 32-bit floats
# Unsure if this is doing too much.
function buggy_generate_gaussians(
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


    #

    
    idx = 1
    for cluster_id in 1:num_centers
        # TODO: Setup standard normal distributions for this cluster_id w/ mean_range and sigma_range
        
        # First randimoly generate Mean for cluster
        cluster_mean = (rand(g) * (mean_range[2] - mean_range[1] +1)) + mean_range[1]
        # Next randimly genearte sigma for cluster
        cluster_sigma = (rand(g) * (sigma_range[2] - sigma_range[1]+1)) + sigma_range[1]
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


function generate_gaussians(
    ::Type{T};
    dim::Int,
    num_centers::Int,
    points_per_center::Int,
    mean_range::Tuple{T, T},
    sigma_range::Tuple{T, T},
    rng::AbstractRNG = Xoshiro(123)
) where {T<:AbstractFloat}
 
    total_points = num_centers * points_per_center
    points = Vector{Vector{T}}(undef, total_points)
    labels = Vector{Int}(undef, total_points)
 
    mean_lo, mean_hi = mean_range
    sigma_lo, sigma_hi = sigma_range
 
    idx = 1
    for cluster_id in 1:num_centers
        # Independent mean per dimension -> center is a random point in R^dim
        cluster_mean = mean_lo .+ rand(rng, dim) .* (mean_hi - mean_lo)
 
        # Independent sigma per dimension -> diagonal covariance
        cluster_sigma = sigma_lo .+ rand(rng, dim) .* (sigma_hi - sigma_lo)
        cov = Diagonal(cluster_sigma .^ 2)
 
        d = MvNormal(cluster_mean, cov)
 
        for p in 1:points_per_center
            points[idx] = rand(rng, d)
            labels[idx] = cluster_id
            idx += 1
        end
    end
 
    return points, labels
end

# figure out how this whole code works
function generate_rings(
    ::Type{T};
    dim::Int,
    num_rings::Int,
    points_per_ring::Int,
    radius_range::Tuple{T, T},
    noise_sigma::T,
    rng::AbstractRNG = Xoshiro(123)
) where {T<:AbstractFloat}
    # 1 -10, 1 - 91 
    total_points = num_rings * points_per_ring
    points = Vector{Vector{T}}(undef, total_points)
    labels = Vector{Int}(undef, total_points)

    radii = num_rings == 1 ?
        T[radius_range[1]] :
        collect(range(radius_range[1], radius_range[2]; length = num_rings)) # figure out what this does

    idx = 1
    for ring_id in 1:num_rings
        r = radii[ring_id]
        
        for p in 1:points_per_ring
            # Uniform random direction on the unit (dim-1)-sphere: draw a
            # standard normal vector and normalize it. This is the standard
            # trick for uniform sampling on a hypersphere surface.
            direction = randn(rng, T, dim)
            direction ./= norm(direction) 

            # Radial noise gives the shell some thickness. Floor at a small
            # positive value so heavy noise can't fold points through the
            # origin and invert which ring they visually belong to.
            noisy_radius = max(r + noise_sigma * randn(rng, T), T(0.01))

            points[idx] = direction .* noisy_radius
            labels[idx] = ring_id
            idx += 1
        end
    end

    return points, labels
end
