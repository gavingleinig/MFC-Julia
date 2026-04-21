# src/Algorithms/MetricForests.jl

struct MFCResult{T}
    cluster_edges::Vector{Tuple{Int, Int, T}}     # tuples representing the edges that form the MST within each local cluster
    completion_edges::Vector{Tuple{Int, Int, T}}   # tuples representing the edges that bridge the separate clusters
    sub_cluster_runtime::Float64                # time (in seconds) taken to compute the local cluster MSTs
    completion_edges_runtime::Float64           # time (in seconds) taken to find the connecting edges
    completion_runtime::Float64                 # total runtime of the MFC algorithm
end

function metric_forest_completion_approx(
    points::AbstractVector{V},          # An array of type V representing data points
    cluster_count::Int,                 # The total number of clusters the data was divided into initial clustering
    cluster_assignments::Vector{Int},   # An array mapping each point to its assigned cluster ID. The length must exactly match the length of 'points.' Values should range from '1' to 'cluster_count'.
    dist_func::F
) where {V, F<:Function}
    
    total_start_time = time()
    
    T = Base.promote_op(dist_func, eltype(points), eltype(points))
    cluster_edges = Tuple{Int, Int, T}[]
    completion_edges = Tuple{Int, Int, T}[]
    
    ### MST Each Cluster ###
    sub_cluster_start = time()
    

    # Ie, points = [X1, X2, X3, X4, X5]
    # Ie, cluster_assignments = [1, 1, 2, 2, 1]

    # Group original point indices by their assigned cluster
    global_indices_by_cluster = [Int[] for _ in 1:cluster_count] # init vectors holding clusters
    for (global_index, cluster_id) in enumerate(cluster_assignments)
        push!(global_indices_by_cluster[cluster_id], global_index)
    end
    # Ie, global_indices_by_cluster = [[1, 2, 5], [3, 4]]
    
    # Run MST on each cluster
    for c in 1:cluster_count
        # Ie, current_cluster_indices = [1, 2, 5]
        current_cluster_indices = global_indices_by_cluster[c]
        
        # A cluster must have at least 2 points to form an edge
        length(current_cluster_indices) < 2 && continue
        
        # Extract data for points in cluster
        # Ie, current_cluster_subset = [X1, X2, X5]
        current_cluster_subset = points[current_cluster_indices]
        
        local_mst_edges = mst_implicit(current_cluster_subset, dist_func)
        
        # Map the local MST indices back to the global dataset indices
        for edge in local_mst_edges
            global_a = current_cluster_indices[edge.a]
            global_b = current_cluster_indices[edge.b]
            edge_weight = edge.weight
            push!(cluster_edges, (global_a, global_b, edge_weight))
        end
    end
    
    sub_cluster_runtime = time() - sub_cluster_start

    ### Completeion Edges ###
    completion_edges_start = time()
    
    # TODO: Implement the completion phase.
    # Ie, find shortest edges between the disconnected clusters to form global MST, or representatives, etc...
    unmapped_completion_edges = metric_forest_completion_Edges_approx_simple(points,cluster_count,global_indices_by_cluster, dist_func)
    unmapped_completion_edges = mst_complete!(cluster_count, unmapped_completion_edges)

    completion_edges_runtime = time() - completion_edges_start
    completion_edges = convert_completion_to_weight(unmapped_completion_edges)

    
    total_runtime = time() - total_start_time

    return MFCResult(
        cluster_edges, 
        completion_edges, 
        sub_cluster_runtime, 
        completion_edges_runtime, 
        total_runtime
    ) 
end

function metric_forest_completion_optimal(
    points::AbstractVector{V},          # An array of type V representing data points
    cluster_count::Int,                 # The total number of clusters the data was divided into initial clustering
    cluster_assignments::Vector{Int},   # An array mapping each point to its assigned cluster ID. The length must exactly match the length of 'points.' Values should range from '1' to 'cluster_count'.
    dist_func::F
) where {V, F<:Function}
    
    total_start_time = time()
    
    T = Base.promote_op(dist_func, eltype(points), eltype(points))
    cluster_edges = Tuple{Int, Int, T}[]
    completion_edges = Tuple{Int, Int, T}[]
    
    ### MST Each Cluster ###
    sub_cluster_start = time()
    

    # Ie, points = [X1, X2, X3, X4, X5]
    # Ie, cluster_assignments = [1, 1, 2, 2, 1]

    # Group original point indices by their assigned cluster
    global_indices_by_cluster = [Int[] for _ in 1:cluster_count] # init vectors holding clusters
    for (global_index, cluster_id) in enumerate(cluster_assignments)
        push!(global_indices_by_cluster[cluster_id], global_index)
    end
    # Ie, global_indices_by_cluster = [[1, 2, 5], [3, 4]]
    
    # Run MST on each cluster
    for c in 1:cluster_count
        # Ie, current_cluster_indices = [1, 2, 5]
        current_cluster_indices = global_indices_by_cluster[c]
        
        # A cluster must have at least 2 points to form an edge
        length(current_cluster_indices) < 2 && continue
        
        # Extract data for points in cluster
        # Ie, current_cluster_subset = [X1, X2, X5]
        current_cluster_subset = points[current_cluster_indices]
        
        local_mst_edges = mst_implicit(current_cluster_subset, dist_func)
        
        # Map the local MST indices back to the global dataset indices
        for edge in local_mst_edges
            global_a = current_cluster_indices[edge.a]
            global_b = current_cluster_indices[edge.b]
            edge_weight = edge.weight
            push!(cluster_edges, (global_a, global_b, edge_weight))
        end
    end
    
    sub_cluster_runtime = time() - sub_cluster_start

    ### Completeion Edges ###
    completion_edges_start = time()
    
    # TODO: Implement the completion phase.
    # Ie, find shortest edges between the disconnected clusters to form global MST, or representatives, etc...
    unmapped_completion_edges = metric_forest_completion_Edges_approx_simple(points,cluster_count,global_indices_by_cluster, dist_func)
    unmapped_completion_edges = mst_complete!(cluster_count, unmapped_completion_edges)

    completion_edges_runtime = time() - completion_edges_start
    completion_edges = metric_forest_completion_edges_optimal(unmapped_completion_edges)

    
    total_runtime = time() - total_start_time

    return MFCResult(
        cluster_edges, 
        completion_edges, 
        sub_cluster_runtime, 
        completion_edges_runtime, 
        total_runtime
    ) 
end


function metric_forest_completion_simple(
    points::AbstractVector{V},          # An array of type V representing data points
    cluster_count::Int,                 # The total number of clusters the data was divided into initial clustering
    cluster_assignments::Vector{Int},   # An array mapping each point to its assigned cluster ID. The length must exactly match the length of 'points.' Values should range from '1' to 'cluster_count'.
    dist_func::F
) where {V, F<:Function}
    
    total_start_time = time()
    
    T = Base.promote_op(dist_func, eltype(points), eltype(points))
    cluster_edges = Tuple{Int, Int, T}[]
    completion_edges = Tuple{Int, Int, T}[]
    
    ### MST Each Cluster ###
    sub_cluster_start = time()
    

    # Ie, points = [X1, X2, X3, X4, X5]
    # Ie, cluster_assignments = [1, 1, 2, 2, 1]

    # Group original point indices by their assigned cluster
    global_indices_by_cluster = [Int[] for _ in 1:cluster_count] # init vectors holding clusters
    for (global_index, cluster_id) in enumerate(cluster_assignments)
        push!(global_indices_by_cluster[cluster_id], global_index)
    end
    # Ie, global_indices_by_cluster = [[1, 2, 5], [3, 4]]
    
    # Run MST on each cluster
    for c in 1:cluster_count
        # Ie, current_cluster_indices = [1, 2, 5]
        current_cluster_indices = global_indices_by_cluster[c]
        
        # A cluster must have at least 2 points to form an edge
        length(current_cluster_indices) < 2 && continue
        
        # Extract data for points in cluster
        # Ie, current_cluster_subset = [X1, X2, X5]
        current_cluster_subset = points[current_cluster_indices]
        
        local_mst_edges = mst_implicit(current_cluster_subset, dist_func)
        
        # Map the local MST indices back to the global dataset indices
        for edge in local_mst_edges
            global_a = current_cluster_indices[edge.a]
            global_b = current_cluster_indices[edge.b]
            edge_weight = edge.weight
            push!(cluster_edges, (global_a, global_b, edge_weight))
        end
    end
    
    sub_cluster_runtime = time() - sub_cluster_start

    ### Completeion Edges ###
    completion_edges_start = time()
    
    # TODO: Implement the completion phase.
    # Ie, find shortest edges between the disconnected clusters to form global MST, or representatives, etc...
    unmapped_completion_edges = metric_forest_completion_edges_random_connection(points,cluster_count,global_indices_by_cluster, dist_func)
    unmapped_completion_edges = mst_complete!(cluster_count, unmapped_completion_edges)

    completion_edges_runtime = time() - completion_edges_start
    completion_edges = metric_forest_completion_edges_optimal(unmapped_completion_edges)

    
    total_runtime = time() - total_start_time

    return MFCResult(
        cluster_edges, 
        completion_edges, 
        sub_cluster_runtime, 
        completion_edges_runtime, 
        total_runtime
    ) 
end