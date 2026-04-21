


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

