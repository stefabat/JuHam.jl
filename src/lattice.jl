
abstract type Topology end

"""
    Chain(L, a, D)

1-dimensional chain with equally spaced sites

# Arguments
* `L::Int`          : number of sites
* `a::AbstractFloat`: lattice constant
* `D::Matrix{Int}`  : min path connectivty between each pair of sites
"""
struct Chain <: Topology
    L::Int
    a::AbstractFloat
    D::Matrix{Int}

    # Constructor
    function Chain(L::Int, a::Real)
        assert(L > 1)
        D = fill(L,(L,L)) 
        D[1:L+1:end]   = 0
        D[2:L+1:end]   = 1
        D[L+1:L+1:end] = 1
        floydwarshall!(D)
        new(L, a, D)
    end

end


"""
    Lattice{N}(L, a, D)

regular `N`-dimensional lattice

# Arguments
* `L::NTuple{N,Integer}`      : number of sites for each of the `N` dimensions
* `a::NTuple{N,AbstractFloat}`: lattice constant in each of the `N` dimensions
* `D::Matrix{Integer}`        : min path connectivty between each pair of sites
"""
struct Lattice{N} <: Topology
    L::NTuple{N,Integer}
    a::NTuple{N,AbstractFloat}
    D::Matrix{Integer}
end
