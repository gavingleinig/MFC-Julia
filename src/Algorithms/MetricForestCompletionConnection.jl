

struct completion_edge
    a::Int
    b::Int
    a_rep::Int
    b_rep::Int
    weight::Float64
end

function metric_forest_completion_Edges_approx_simple(
    points::AbstractVector{V},          # An array of type V representing data points
    cluster_count::Int,                 # The total number of clusters the data was divided into initial clustering
    global_indices_by_cluster::Vector{Vector{Int}},   # An array mapping each point to its assigned cluster ID. The length must exactly match the length of 'points.' Values should range from '1' to 'cluster_count'.

    dist_func::F
) where {V, F<:Function}
    unmapped_completion_edges = Vector{completion_edge}()
    size_hint!(unmapped_completion_edges, cluster_count * (cluster_cout - 1))

    # Consider all clusters 
    for i in 1:cluster_count
        partition_i = global_indices_by_cluster[i]
        for j in 2:cluster_count
            partition_j = global_indices_by_cluster[j]

            if length(partition_i) == 0 || length(partition_j) ==0
                continue
            end

            # For simple we wil use 0 for our repersentives
            partition_i_rep = 0
            partition_j_rep = 0

            current_completion_edges = completion_edge(i,j,0,0,Inf)
            

            for partition_j_edge_num in partition_j
                dist = dist_func(points[global_indices_by_cluster[i][i_rep]], points[partition_j_edge_num])
                if dist < current_completion_edges.weight
                    current_completion_edges.a_rep = global_indices_by_cluster[i][i_rep]
                    current_completion_edges.b_rep = partition_j_edge_num
                    current_completion_edges.weight = dist
                end
            end

            for partition_i_edge_num in partition_i
                dist = dist_func(points[partition_i_edge_num],points[global_indices_by_cluster[j][j_rep]])
                if dist < current_completion_edges.weight
                    current_completion_edges.a_rep = global_indices_by_cluster[i][i_rep]
                    current_completion_edges.b_rep = partition_j_edge_num
                    current_completion_edges.weight = dist
                end
            end

            push!(unmapped_completion_edges, current_completion_edges)

        end

        return unmapped_completion_edges
    end





    # Choose





end