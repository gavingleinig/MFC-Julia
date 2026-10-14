
using DataStructures
struct CompletionEdge{T<:Real} 
    a::Int
    b::Int
    a_rep::Int
    b_rep::Int
    weight::T
end

struct RepInfo
    index::Int
    distance::Float64
end


function metric_forest_completion_Edges_approx_simple(
    points::AbstractVector{V},          # An array of type V representing data points
    cluster_count::Int,                 # The total number of clusters the data was divided into initial clustering
    global_indices_by_cluster::Vector{Vector{Int}},  

    dist_func::F
) where {V, F<:Function}
    T = Base.promote_op(dist_func, eltype(points), eltype(points))
    unmapped_completion_edges = Vector{CompletionEdge{T}}()
    # size_hint!(unmapped_completion_edges, cluster_count * (cluster_cout - 1))

    # Consider all clusters 
    for i in 1:cluster_count
        partition_i = global_indices_by_cluster[i]

        for j in i+1:cluster_count
            partition_j = global_indices_by_cluster[j]

            if length(partition_i) == 0 || length(partition_j) == 0
                continue
            end

            # For simple we wil use first node in cluster for our repersentives
            i_rep = 1
            j_rep = 1

            best_dist = Inf
            best_a_rep = -1
            best_b_rep = -1

            for point_j_idx in partition_j
                dist = dist_func(points[global_indices_by_cluster[i][i_rep]], points[point_j_idx])
                    # Maybe move points[global_indices_by_cluster[i][i_rep]] up out of this loop
                    # So like:
                    # rep_idx_i = partition_i[1] 
                    # point_i_rep = points[rep_idx_i]
                    # 
                    # because: point_i_rep == points[global_indices_by_cluster[i][i_rep]]?
                    # So then: dist = dist_func(point_i_rep, points[point_j_idx])
                if dist < best_dist
                    best_dist = dist
                    best_a_rep = global_indices_by_cluster[i][i_rep] # = rep_idx_i
                    best_b_rep = point_j_idx
                end
            end

            for point_i_idx in partition_i
                dist = dist_func(points[point_i_idx],points[global_indices_by_cluster[j][j_rep]])
                # Similarly,
                # rep_idx_j = partition_j[1]
                # point_j_rep = points[rep_idx_j]
                # dist = dist_func(points[point_i_idx], point_j_rep)
                if dist < best_dist
                    best_dist = dist
                    best_a_rep = point_i_idx
                    best_b_rep = global_indices_by_cluster[j][j_rep] # = rep_idx_j
                end
            end

            push!(
                unmapped_completion_edges, 
                CompletionEdge(i, j, best_a_rep, best_b_rep, best_dist)
            )
        end 
    end

    return unmapped_completion_edges
    # Choose
end


function metric_forest_completion_edges_optimal(
    points::AbstractVector{V},          # An array of type V representing data points
    cluster_count::Int,                 # The total number of clusters the data was divided into initial clustering
    global_indices_by_cluster::Vector{Vector{Int}},  

    dist_func::F
) where {V, F<:Function}
    T = Base.promote_op(dist_func, eltype(points), eltype(points))
    unmapped_completion_edges = Vector{CompletionEdge{T}}()

    # Consider all clusters 
    for cluster_i in 1:cluster_count
        partition_i = global_indices_by_cluster[cluster_i]

        for cluster_j in (cluster_i+1):cluster_count
            partition_j = global_indices_by_cluster[cluster_j]

            if length(partition_i) == 0 || length(partition_j) == 0
                continue
            end

            # For simple we wil use first node in cluster for our repersentives
            i_rep = 1
            j_rep = 1

            best_dist = Inf
            best_a_rep = -1
            best_b_rep = -1

            for i in 1:length(global_indices_by_cluster[cluster_i])
                for j in 1:length(global_indices_by_cluster[cluster_j])
                    a_rep = global_indices_by_cluster[cluster_i][i]
                    b_rep = global_indices_by_cluster[cluster_j][j]
                    dist = dist_func(points[a_rep],points[b_rep])
                    if best_dist > dist
                        best_dist = dist

                        best_a_rep = a_rep
                        best_b_rep = b_rep
                    end
                end
            end

           

            push!(
                unmapped_completion_edges, 
                CompletionEdge(cluster_i, cluster_j, best_a_rep, best_b_rep, best_dist)
            )
        end 
    end

    return unmapped_completion_edges


end

function metric_forest_completion_edges_random_connection(
    points::AbstractVector{V},          # An array of type V representing data points
    cluster_count::Int,                 # The total number of clusters the data was divided into initial clustering
    global_indices_by_cluster::Vector{Vector{Int}},  

    dist_func::F
) where {V, F<:Function}
    T = Base.promote_op(dist_func, eltype(points), eltype(points))
    unmapped_completion_edges = Vector{CompletionEdge{T}}()

    # Consider all clusters 
    for i in 1:cluster_count
        partition_i = global_indices_by_cluster[i]

        for j in (i+1):cluster_count
            partition_j = global_indices_by_cluster[j]

            if length(partition_i) == 0 || length(partition_j) == 0
                continue
            end

            # For simple we wil use first node in cluster for our repersentives
            i_rep = 1
            j_rep = 1

            dist = dist_func(points[global_indices_by_cluster[i][i_rep]],points[global_indices_by_cluster[j][j_rep]])

            

            push!(
                unmapped_completion_edges, 
                CompletionEdge(i, j, global_indices_by_cluster[i][i_rep], global_indices_by_cluster[j][j_rep], dist)
            )
        end 
    end

    return unmapped_completion_edges


end



# Multi repersintative completion edges

function get_unmapped_completion_edges_from_reps(
    points::AbstractVector{V},
    cluster_vecs::Vector{Vector{Int}},
    rep_vecs::Vector{Vector{Tuple{Int,Float32}}},
    dist_func::F
) where {V, F<:Function}

    cluster_count = length(cluster_vecs)
    T = Base.promote_op(dist_func, eltype(points), eltype(points))

    unmapped_completion_edges = Vector{CompletionEdge{T}}()
    sizehint!(unmapped_completion_edges, cluster_count * (cluster_count - 1))

    runtime = @elapsed begin
        for i in 1:(cluster_count - 1)
            for j in 2:cluster_count

                if isempty(cluster_vecs[i]) || isempty(cluster_vecs[j])
                    continue
                end

                best_weight = Inf
                best_a_rep = -1
                best_b_rep = -1

                # every rep of cluster i vs every point of cluster j
                for point_j_idx in cluster_vecs[j]
                    for (i_rep, _cost) in rep_vecs[i]
                        dist = dist_func(points[i_rep], points[point_j_idx])
                        if dist < best_weight
                            best_weight = dist
                            best_a_rep = i_rep
                            best_b_rep = point_j_idx
                        end
                    end
                end

                # every point of cluster i vs every rep of cluster j
                for point_i_idx in cluster_vecs[i]
                    for (j_rep, _cost) in rep_vecs[j]
                        dist = dist_func(points[point_i_idx], points[j_rep])
                        if dist < best_weight
                            best_weight = dist
                            best_a_rep = point_i_idx
                            best_b_rep = j_rep
                        end
                    end
                end

                push!(
                    unmapped_completion_edges,
                    CompletionEdge(i, j, best_a_rep, best_b_rep, best_weight)
                )
            end
        end
    end

    return unmapped_completion_edges, runtime
end
