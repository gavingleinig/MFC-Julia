# src/Algorithms/MetricForests.jl

struct MFCResult
    cluster_edges::Vector{Tuple{Int, Int}}      # tuples representing the edges that form the MST within each local cluster
    completion_edges::Vector{Tuple{Int, Int}}   # tuples representing the edges that bridge the separate clusters
    sub_cluster_runtime::Float64                # time (in seconds) taken to compute the local cluster MSTs
    completion_edges_runtime::Float64           # time (in seconds) taken to find the connecting edges
    completion_runtime::Float64                 # total runtime of the MFC algorithm
end

function metric_forest_completion(
    points::AbstractVector{V},          # An array of type V representing data points
    cluster_count::Int,                 # The total number of clusters the data was divided into initial clustering
    cluster_assignments::Vector{Int},   # An array mapping each point to its assigned cluster ID. The length must exactly match the length of 'points.' Values should range from '1' to 'cluster_count'.
    dist_func::F
) where {V, F<:Function}
    
    total_start_time = time()
    
    cluster_edges = Tuple{Int, Int}[]
    completion_edges = Tuple{Int, Int}[]
    
    ### MST Each Cluster ###
    sub_cluster_start = time()
    
    # Group original point indices by their assigned cluster
    cluster_indices = [Int[] for _ in 1:cluster_count] # init vectors holding clusters
    for (global_index, cluster_id) in enumerate(cluster_assignments)
        push!(cluster_indices[cluster_id], global_index)
    end
    
    # Run MST on each cluster
    for c in 1:cluster_count
        local_indices = cluster_indices[c]
        
        # A cluster must have at least 2 points to form an edge
        length(local_indices) < 2 && continue
        
        # Extract data for points in cluster
        cluster_points = points[local_indices]
        
        local_mst_edges = mst_implicit(cluster_points, dist_func)
        
        # Map the local MST indices back to the global dataset indices
        for edge in local_mst_edges
            global_a = local_indices[edge.a]
            global_b = local_indices[edge.b]
            push!(cluster_edges, (global_a, global_b))
        end
    end
    
    sub_cluster_runtime = time() - sub_cluster_start

    ### Completeion Edges ###
    completion_edges_start = time()
    
    # TODO: Implement the completion phase.
    # Ie, find shortest edges between the disconnected clusters to form global MST, or representatives, etc...
    
    completion_edges_runtime = time() - completion_edges_start
    
    total_runtime = time() - total_start_time

    return MFCResult(
        cluster_edges, 
        completion_edges, 
        sub_cluster_runtime, 
        completion_edges_runtime, 
        total_runtime
    ) 
end