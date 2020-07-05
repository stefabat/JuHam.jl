
using LinearAlgebra

"""
    floydwarshall(A::Matrix{Int})

Floyd-Warshall algorithm to find all-pairs shortest paths in O(n^3) time.

`A` is the NxN adjacency matrix of the graph, with zeros on the diagonal,
one for every pair of connecting vertices and at least `N` for the rest.
"""
function floydwarshall!(A::Matrix{Int})
    @assert(issymmetric(A))
    N = size(A,1)
    for k = 1:N
        for i = 1:N
            for j = 1:N
                A[i,j] = min(A[i,j],A[i,k] + A[k,j])
            end
        end
    end
    return A
end