using DataStructures

struct WeightedEdge{T<:Real} # T must be a subtype of Real, ie Int, Float64, etc
    weight::T
    a::Int
    b::Int
end



# Tell Julia how to compare WeigthedEdges; so can sort() later on
Base.isless(e1::WeightedEdge, e2::WeightedEdge) = isless(e1.weight, e2.weight)
Base.isless(e1::CompletionEdge, e2::CompletionEdge) = isless(e1.weight, e2.weight)
# Compute Optimum MST for given graph
function mst!(n::Int, edges::Vector{WeightedEdge{T}}) where T
    sort!(edges) # uses above comparison, Base.isless
    dj_sets = IntDisjointSets(n)
    kept_edges = WeightedEdge{T}[]
    
    for e in edges
        if !in_same_set(dj_sets, e.a, e.b) # If edge connects two verticies in different trees:
            union!(dj_sets, e.a, e.b) # merge two trees
            push!(kept_edges, e) # add edge to MST
            
            # Stop when there are n-1 edges in MST
            length(kept_edges) == n - 1 && break
        end 
    end
    return kept_edges
end

# Takes list of points and distance function, returns list of edges with distances pre-calculated
function mst_implicit(points::AbstractVector{V}, dist_func::F) where {V, F}
    n = length(points)
    T = Base.promote_op(dist_func, eltype(points), eltype(points)) # Infer return type of distance function

    if n < 2
        return WeightedEdge{T}[] 
    end
    
    edges = WeightedEdge{T}[]
    
    # Batched implicit generation
    for i in 1:(n-1)
        for j in (i+1):n
            w = dist_func(points[i], points[j])
            push!(edges, WeightedEdge(w, i, j))
        end
        
        
        if length(edges) > 10_000_000 # Batch threshold 
            # Compute MST and prune
            edges = mst!(n, edges)
        end
    end
    
    # Final MST calculation
    return mst!(n, edges)
end


function mst_complete!(n::Int, edges::Vector{CompletionEdge{T}}) where T
    sort!(edges) # uses above comparison, Base.isless
    dj_sets = IntDisjointSets(n)
    kept_edges = CompletionEdge{T}[]
    
    for e in edges
        if !in_same_set(dj_sets, e.a, e.b) # If edge connects two verticies in different trees:
            union!(dj_sets, e.a, e.b) # merge two trees
            push!(kept_edges, e) # add edge to MST
            
            # Stop when there are n-1 edges in MST
            length(kept_edges) == n - 1 && break
        end 
    end
    return kept_edges
end
