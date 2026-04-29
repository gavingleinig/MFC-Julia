


function convert_completion_to_weight(unmapped_completion_edges::Vector{CompletionEdge{T}}) where T
    completion_edges = Vector{Tuple{Int, Int, T}}()
    for e in unmapped_completion_edges
        push!(completion_edges, (e.a_rep,e.b_rep,e.weight))
    end



    return completion_edges

end


function convert_mst_to_weight(MST_edges::Vector{WeightedEdge{T}}) where T
    completion_edges = Vector{Tuple{Int, Int, T}}()
    for e in MST_edges
        push!(completion_edges, (e.a,e.b,e.weight))
    end



    return completion_edges

end


"""
Calculates the practical upper bound for the y-overlap (gamma).
This measures how much the optimal MST overlaps with the initial k-center partitions.
"""
function compute_gamma_overlap(
    initial_forest_edges::Vector{Tuple{Int, Int, T}}, 
    optimal_mst_edges::Vector{Tuple{Int, Int, T}},
    assignments::Vector{Int}
) where T <: Real
    # Weight of the initial forest
    w_initial_forest = sum(e[3] for e in initial_forest_edges)

    # Weight of the optimal MST  within the initial clusters
    w_optimal_internal = zero(T)
    for (u, v, w) in optimal_mst_edges
        if assignments[u] == assignments[v]
            w_optimal_internal += w
        end
    end

    if w_optimal_internal == zero(T)
        return Inf
    end

    return w_initial_forest / w_optimal_internal
end