using Distances
using Graphs

println("Hello, World!")

# Test Distances.jl
p1 = [1.0, 2.0]
p2 = [4.0, 6.0]
dist = evaluate(Euclidean(), p1, p2)
println("Euclidean distance between $p1 and $p2 is: $dist")

# Test Graphs.jl
g = path_graph(4)
println("Created a simple graph with $(nv(g)) vertices and $(ne(g)) edges.")