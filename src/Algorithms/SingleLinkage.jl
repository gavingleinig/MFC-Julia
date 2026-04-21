using Graphs
using Clustering
# Could also do this w/o Graphs package using a Disjoint Set data structure?

# Helper fctn to convert a list of connected components into a flat assignments array
function components_to_assignments(components::Vector{Vector{Int}}, num_points::Int)
    assignments = Vector{Int}(undef, num_points) # allocate memory, don't init values
    for (cluster_id, comp) in enumerate(components)
        for node in comp
            assignments[node] = cluster_id
        end
    end
    return assignments
end

function single_linkage_threshold(
    edges::Vector{Tuple{Int, Int, T}}, 
    num_points::Int, 
    threshold::T
) where T<:Real
    start_time = time()
    
    # Create an empty graph with num_points vertices and 0 edges
    g = SimpleGraph(num_points)
    
    # Add only the edges within threshold
    for (u, v, w) in edges
        if w <= threshold
            add_edge!(g, u, v)
        end
    end
    
    # Find the isolated groups, our clusters
    comps = connected_components(g)
    assignments = components_to_assignments(comps, num_points)
    
    return ClusteringResult{Float64}(assignments, time() - start_time)
end

function single_linkage_k_clusters(
    edges::Vector{Tuple{Int, Int, T}}, 
    num_points::Int, 
    k::Int
) where T<:Real
    start_time = time()
    
    # Sort edges by weight, ascending
    sorted_edges = sort(edges, by = x -> x[3])
    
    # A ST with N vertices and N - k edges has k connected components
    edges_to_keep = num_points - k
    
    # Make suer not keeping negative edges, ie k > num_points
    edges_to_keep = max(0, edges_to_keep) 
    
    g = SimpleGraph(num_points)
    
    # Add only the smallest edges_to_keep edges
    for i in 1:min(edges_to_keep, length(sorted_edges))
        u, v, _ = sorted_edges[i]
        add_edge!(g, u, v)
    end
    
    comps = connected_components(g)
    assignments = components_to_assignments(comps, num_points)
    
    return ClusteringResult{Float64}(assignments, time() - start_time)
end

function best_single_linkage_threshold(
    edges::Vector{Tuple{Int, Int, Float64}}, 
    num_points::Int,
    ground_truth::AbstractVector{Int64},
    increment::Float64
)
    max = edges[1]
    min = edges[1]

    for e in edges
        weight = e[3]

        if max[3] < weight
            max = e
        end

        if min[3] > weight
            min = e
        end
    end

    
    current = min[3]


    ari_max = 0
    best_threshold = single_linkage_threshold(edges, num_points, current)


    while current < max[3]
        threshold_result = single_linkage_threshold(edges, num_points,current)


        ari = randindex(threshold_result.assignments, ground_truth)[1]
        if ari_max < ari
            ari_max = ari
            best_threshold = threshold_result
        end

        


        current += increment
    end

    return (best_threshold, ari_max)
end