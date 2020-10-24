
abstract type Topology end


"""
    Chain <: Topology

Type representing a 1-dimensional chain of length `L` with equally spaced sites
and lattice constant `a`.
"""
struct Chain <: Topology
    L::Int              # number of sites
    a::Float64          # lattice constant (intersite spacing)
    D::Matrix{Int}      # connectivity matrix
end


"""
    Chain(L::Int, a::Real)

Construct a 1-dimensional `Chain` of length `L` and intersite distance `a`.

`L` must be an integer number greater than one. The intersite distance `a`,
that is, the lattice constant, must be a positive real number.
"""
function Chain(L::Int, a::Real)
    if L <= 1
        error("`L` should be an integer number greater than one.")
    end
    if a <= 0
        error("`a` should be a real number greater than zero.")
    end
    # construct the connectivity matrix
    D = fill(L,(L,L)) 
    D[1:L+1:end]   .= 0
    D[2:L+1:end]   .= 1
    D[L+1:L+1:end] .= 1
    # determine the connectivity
    floydwarshall!(D)
    
    return Chain(L, a, D)
end


# """
#     Lattice{N}(L, a, D)

# regular `N`-dimensional lattice

# # Arguments
# * `L::NTuple{N,Integer}`      : number of sites for each of the `N` dimensions
# * `a::NTuple{N,AbstractFloat}`: lattice constant in each of the `N` dimensions
# * `D::Matrix{Integer}`        : min path connectivty between each pair of sites
# """
# struct Lattice{N} <: Topology
#     L::NTuple{N,Integer}
#     a::NTuple{N,AbstractFloat}
#     D::Matrix{Integer}
# end


"""
    floydwarshall(A::Matrix{Int})

Floyd-Warshall algorithm to find the shortest path between each node of a graph in O(n^3) time.

`A` is the NxN adjacency matrix of the graph, with zeros on the diagonal,
one for every pair of connecting vertices and a value equal or greater than
N for all remaining elements.
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
