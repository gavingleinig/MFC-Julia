using DataStructures

struct WeightedEdge{T<:Real} # T must be a subtype of Real, ie Int, Float64, etc
    weight::T
    a::Int
    b::Int
end

# Tell Julia how to compare WeigthedEdges; so can sort() later on
Base.isless(e1::WeightedEdge, e2::WeightedEdge) = isless(e1.weight, e2.weight)
Base.isless(e1::CompletionEdge, e2::CompletionEdge) = isless(e1.weight, e2.weight)


# Computes optimum MST for given graph using Kruskal's algorithm
function mst!(n::Int, edges::Vector{E}) where E
    sort!(edges) # uses above comparison, Base.isless
    dj_sets = IntDisjointSets(n)
    kept_edges = E[]
    
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
    beginTime = time()
    n = length(points)
    T = Base.promote_op(dist_func, eltype(points), eltype(points)) # Infer return type of distance function

    if n < 2
        return WeightedEdge{T}[] 
    end
    
    edges = WeightedEdge{T}[]
    
    # Batched implicit generation of all edges

    for i in 1:(n-1) # for every point
        for j in (i+1):n # compare with every other points
            w = dist_func(points[i], points[j])
            push!(edges, WeightedEdge(w, i, j))
        end
        
        if length(edges) > 50_000_000 # Batch threshold 
            # Compute MST and prune
            edges = mst!(n, edges)
        end
    end
    
    # Final MST calculation


    final_edges = mst!(n, edges)


    runtime = time() - beginTime

    return final_edges
end

# Computes a naive, random spanning tree by randomly selecting edges 
# and adding them if they do not form a cycle.
function naive_random_st(points::AbstractVector{V}, dist_func::F) where {V, F}
    n = length(points)
    T = Base.promote_op(dist_func, eltype(points), eltype(points)) # Infer return type

    if n < 2
        return WeightedEdge{T}[] 
    end

    dj_sets = IntDisjointSets(n)
    kept_edges = WeightedEdge{T}[]
    sizehint!(kept_edges, n - 1) # Optimize memory allocation
    
    # Continually pick random pairs until we have n-1 valid edges
    while length(kept_edges) < n - 1
        a = rand(1:n)
        b = rand(1:n)
        
        # If the vertices are distinct and belong to different trees
        if a != b && !in_same_set(dj_sets, a, b)
            union!(dj_sets, a, b) # Merge the two trees
            
            # Compute distance only when we know we are keeping the edge
            w = dist_func(points[a], points[b])
            push!(kept_edges, WeightedEdge(w, a, b))
        end
    end
    
    return kept_edges
end