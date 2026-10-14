
struct RepInfo
    index::Int
    distance::Float64
end


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

function convert_mst_to_weight(edges::Vector{CompletionEdge{T}}) where T
    weighted = [WeightedEdge(e.weight, e.a_rep, e.b_rep) for e in edges]
    return convert_mst_to_weight(weighted)
end


"""
    compute_gamma_overlap(initial_forest_edges, optimal_mst_edges, assignments)

Calculate the practical upper bound for the y-overlap. 
This is the weight ratio of the initial local forest to the optimal MST contained within those same clusters
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





# ----------------------------------------------------------------------
# get_best_reps: greedy farthest-point sampling (k-center greedy) to build
# a ranked, cost-monotonic list of representatives for one cluster.
# This is what produces the `reps` vectors that feed `pick_reps` above.
# ----------------------------------------------------------------------

function get_best_reps(
    cluster_count::Integer, 
    points::AbstractVector, 
    cluster::Vector{<:Integer},
    amount::Integer, 
    dist_func::F
)where{V, F<:Function}


    amount <= 0 && return Tuple{Int,Float32}[]

    n = length(cluster)
    amount = min(amount, n)  # can't have more distinct reps than points in the cluster

    # Materialize the cluster's points once: turns `points[cluster[i]]`
    # (two indirections, called amount * n times) into one indirection
    # per call, at a one-time O(n) copy cost.
    local_pts = points[cluster]

    rep_local = Vector{Int}(undef, amount)      # index into `cluster` / `local_pts`
    rep_costs = Vector{Float32}(undef, amount)

    rep_local[1] = 1

    cur_distances = Vector{Float32}(undef, n)
    cur_distances[1] = 0.0f0

    second_local = 1
    max_dist = 0.0f0

    @inbounds for i in 2:n
        d = Float32(dist_func(local_pts[i], local_pts[1]))
        cur_distances[i] = d
        if d > max_dist
            max_dist = d
            second_local = i
        end
    end
    rep_costs[1] = max_dist

    if amount == 1
        return [(cluster[rep_local[1]], rep_costs[1])]
    end

    rep_local[2] = second_local

    # Incrementally add reps 3..amount. Each iteration finalizes the cost
    # of the most-recently-added rep, then finds the next farthest point.
    k = 2
    while true
        last_local = rep_local[k]

        cur_distances[1] = min(cur_distances[1], Float32(dist_func(local_pts[1], local_pts[last_local])))
        max_dist = cur_distances[1]
        new_local = 1

        @inbounds for i in 2:n
            d = Float32(dist_func(local_pts[i], local_pts[last_local]))
            if d < cur_distances[i]
                cur_distances[i] = d
            end
            if cur_distances[i] > max_dist
                max_dist = cur_distances[i]
                new_local = i
            end
        end

        rep_costs[k] = max_dist

        k == amount && break

        k += 1
        rep_local[k] = new_local
    end

    return [(cluster[rep_local[i]], rep_costs[i]) for i in 1:amount]
end

# ----------------------------------------------------------------------
# Stubs for the external helpers referenced in the original snippet.
# Replace these with your real implementations / bindings.
# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
# Core DP: pick how many extra representatives to take from each cluster
# ----------------------------------------------------------------------

"""
    pick_reps(all_reps, cluster_count, budget)

Given `all_reps` (one entry per cluster, each holding a cost-per-rep-count
table sorted cheapest-first) and a shared `budget` of extra representatives,
find the allocation across clusters that minimizes total cost.

Returns `(final_reps, final_cost, work)`:
  - `final_reps[t]`  : chosen reps (mandatory first rep + any extra ones)
                       for cluster `t`
  - `final_cost`     : total cost of the optimal allocation
  - `work`           : (B+1) x T matrix; `work[b+1, t]` = number of reps
                       given to cluster `t` under budget `b`
"""
function pick_reps(all_reps::Vector, cluster_count::Integer, budget::Integer)
    T = cluster_count
    B = budget

    # Clamp any budget-derived index to what this cluster's rep list actually has.
    @inline clamp_idx(reps, i) = min(i + 1, length(reps))

    final_reps = [Vector{eltype(all_reps[t].reps)}() for t in 1:T]

    for t in 1:T
        push!(final_reps[t], all_reps[t].reps[1])
    end

    final_cost = sum(c.reps[1][2] for c in all_reps)

    if B <= 0
        return final_reps, final_cost, nothing
    end

    costs = zeros(Float64, B + 1)
    work  = zeros(Int, B + 1, T)

    reps_1 = all_reps[1].reps
    for b in 0:B
        costs[b + 1]   = reps_1[clamp_idx(reps_1, b)][2]
        work[b + 1, 1] = b
    end

    for t in 2:T
        reps_t = all_reps[t].reps
        new_costs = similar(costs)

        for b in B:-1:0
            best_cost = Inf
            best_i    = 0

            for i in 0:b
                cur_cost = costs[i + 1] + reps_t[clamp_idx(reps_t, b - i)][2]
                if cur_cost < best_cost
                    best_cost = cur_cost
                    best_i    = i
                end
            end

            new_costs[b + 1] = best_cost
            given_to_t       = b - best_i
            work[b + 1, t]   = given_to_t

            for j in 1:(t - 1)
                work[b + 1, j] = work[best_i + 1, j]
            end

            total = sum(@view work[b + 1, 1:t])
            if total != b
                println("t: $t      $total != $b")
            end
        end

        costs .= new_costs
    end

    final_cost = costs[B + 1]

    for t in 1:T
        reps_t = all_reps[t].reps
        for i in 1:work[B + 1, t]
            idx = clamp_idx(reps_t, i)
            idx > length(reps_t) && continue   # nothing further to add for this cluster
            push!(final_reps[t], reps_t[idx])
        end
    end

    return final_reps, final_cost, work
end

# ----------------------------------------------------------------------
# Surrounding logic: timing, sanity check, and building completion edges
# ----------------------------------------------------------------------
