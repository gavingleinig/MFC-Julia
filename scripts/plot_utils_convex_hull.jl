using LinearAlgebra
using Plots

"""
    to_2d(points) -> Vector{Vector{Float64}}

Returns 2D coordinates for plotting. If `points` are already 2D, they're
returned as-is. Otherwise, projects to 2D via PCA (SVD) so higher-dimensional
data (e.g. dim=16 in the sigma sweep) can still be visualized.

Note: convex hulls below are computed on these *projected* 2D coordinates,
not in the original dim-dimensional space. That's an approximation good
enough for visual sanity-checking, not a guarantee that a cluster convex in
full dimension stays convex after projection.
"""
function to_2d(points::Vector{Vector{T}}) where {T<:Real}
    dim = length(points[1])
    if dim == 2
        return [Float64.(p) for p in points]
    end

    X = permutedims(reduce(hcat, points))   # n x dim
    mu = vec(sum(X, dims=1) ./ size(X, 1))
    Xc = X .- mu'
    _, _, V = svd(Xc)
    proj = Xc * V[:, 1:2]                   # n x 2
    return [proj[i, :] for i in 1:size(proj, 1)]
end

"""
    convex_hull(pts) -> Vector{Vector{Float64}}

Andrew's monotone chain convex hull. `pts` is a Vector of 2-element points.
Returns hull vertices in counter-clockwise order. No external dependency
(no need for LazySets/QHull etc).
"""
function convex_hull(pts::AbstractVector{<:AbstractVector{<:Real}})
    n = length(pts)
    n < 3 && return Vector{Float64}[Float64.(p) for p in pts]

    sorted = sort(pts, by = p -> (p[1], p[2]))
    cross(o, a, b) = (a[1] - o[1]) * (b[2] - o[2]) - (a[2] - o[2]) * (b[1] - o[1])

    lower = Vector{Float64}[]
    for p in sorted
        while length(lower) >= 2 && cross(lower[end-1], lower[end], p) <= 0
            pop!(lower)
        end
        push!(lower, Float64.(p))
    end

    upper = Vector{Float64}[]
    for p in reverse(sorted)
        while length(upper) >= 2 && cross(upper[end-1], upper[end], p) <= 0
            pop!(upper)
        end
        push!(upper, Float64.(p))
    end

    return vcat(lower[1:end-1], upper[1:end-1])
end

"""
    plot_clusters_vs_ground_truth(coords2d, ground_truth, predicted; title)

- Fill color of each point  = ground-truth cluster label (from generate_gaussians).
- Black outline polygon     = convex hull of each PREDICTED cluster (from a
  clustering algorithm's assignments), drawn on top.

Comparing the two visually shows where a clustering method's boundaries
diverge from the true generating clusters: same-colored points split across
multiple hulls means the method fragmented a true cluster; a hull containing
multiple colors means it merged separate true clusters together.
"""
function plot_clusters_vs_ground_truth(
    coords2d::Vector{Vector{Float64}},
    ground_truth::Vector{Int},
    predicted::Vector{Int};
    title::String = "Predicted cluster borders vs. ground-truth fill",
    hull_color = :black,
)
    xs = [c[1] for c in coords2d]
    ys = [c[2] for c in coords2d]
    gt_labels = string.("GT ", ground_truth)

    plt = scatter(
        xs, ys; group = gt_labels,
        markersize = 2.5, markerstrokewidth = 0, markeralpha = 0.85,
        legend = false, title = title, titlefontsize = 10,
        xlabel = "dim 1", ylabel = "dim 2",
    )

    for cluster_id in sort(unique(predicted))
        idxs = findall(==(cluster_id), predicted)
        length(idxs) < 3 && continue

        cluster_pts = [coords2d[i] for i in idxs]
        hull = convex_hull(cluster_pts)
        length(hull) < 3 && continue

        hx = [h[1] for h in hull]; push!(hx, hx[1])
        hy = [h[2] for h in hull]; push!(hy, hy[1])
        plot!(plt, hx, hy; seriestype = :path, linecolor = hull_color,
              linewidth = 1.3, linealpha = 0.8, label = "")
    end

    return plt
end