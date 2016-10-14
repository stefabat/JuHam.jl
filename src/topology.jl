## Topology objects

"""
    Topology(dim, L, coords, bonds, dist_mat, bc)

Composite type describing the geometry and the connections of lattices and arbitrary three-dimensional structures.

# Arguments
* `dim::Integer`: Number of active dimensions
* `  L::Integer`: Number of sites/atoms
* `coords::Array{Float64}`: N x dim matrix containing the coordinates of the system
* `bonds::Set{Tuple}`: Set of index pairs indicating directly connected sites/atoms
* `dist_mat::Array{Int64}`: Matrix containing the shortest path length connecting two sites/atoms
* `bc::AbstractString`: String indicating the boundary conditions
"""
type Topology
    dim     ::Integer
    L       ::Integer
    coords  ::Array{Float64}
    bonds   ::Set{Tuple}        # NOTE: deprecated
    dist_mat::Array{Int64,2}
    bc      ::AbstractString    # NOTE: Should be a Hamiltonian property but it's easier to have them here
end


"""
    lattice_generator(L, lat_const = 1.0, bc = "obc")

Generate the topology of a one-dimensional lattice of `L` sites.
"""
function lattice_generator(L::Integer, lat_const::Real = 1.0, bc::String = "obc")
    if bc != "obc" || bc != "pbc"   # Checking that boundary conditions are valid
        throw(DomainError("Possible values are "obc" or "pbc"."))
    end
    coords = collect(linspace(-(L-1)*lat_const/2,(L-1)*lat_const/2,L))  # Generate coordinates of the lattice
    bonds = Set{Tuple}()                # Initialize bonds NOTE: deprecated
    dist_mat = ones(L,L)*(2*L)          # Initialize distance matrix: all entries are 2*L

    # Setting nearest neighbors
    for i = 1:L-1
        push!(bonds,(i,i+1))
        push!(bonds,(i+1,i))
        dist_mat[i,i+1] = 1
        dist_mat[i+1,i] = 1
        dist_mat[i,i]   = 0
    end
    dist_mat[L,L] = 0

    # Add connections if periodic boundary conditions 
    if bc == "pbc"
        dist_mat[1,L] = 1
        dist_mat[L,1] = 1
    end

    # Compute all the distances between the sites
    populate_dist_mat!(dist_mat)

    return Topology(1, L, coords, bonds, dist_mat, bc)
end


"""
    polygon_generator(L, bond_length = 1.0)

Generate the topology of a regular polygon of `L` vertex.

"""
function polygon_generator(L::Integer, bond_length::Real = 1.0)
    if L < 3
        throw(DomainError("L has to be at least 4"))    # Polygons can exist only with at least 3 vertex
    end
    radius = 1 / (2*sin(pi/L))          # Radius of the outer circumference surrounding the polygon
    coords = zeros(L,2)                 # Initialize coords matrix, just 2 dimensions

    # Put the L sites equidistantly on a circle
    coords[:,1] = radius * cos((2*pi/L)*collect(0:(L-1)))
    coords[:,2] = radius * sin((2*pi/L)*collect(0:(L-1)))

    # Setting nearest neighbors
    dist_mat = ones(L,L)*(2*L)          # Initialize distance matrix
    bonds = Set{Tuple}()                # Initialize bonds NOTE: deprecated
    for i = 1:L-1                       # Push bonds into list and set 0 and 1 distances
        push!(bonds,(i,i+1))
        push!(bonds,(i+1,i))
        dist_mat[i,i+1] = 1
        dist_mat[i+1,i] = 1
        dist_mat[i,i]   = 0
    end
    dist_mat[L,L] = 0
    dist_mat[1,L] = 1; dist_mat[L,1] = 1

    # Populate all the distances between the sites
    populate_dist_mat!(dist_mat)

    return Topology(2, L, coords, bonds, dist_mat, "obc")
end


"""
    graphene_generator(L, C, CC_bond = 1.422, bc = "obc")

Generate the topology of a graphene nanoribbon by stacking `C` polyene chains of length `L` on top of each other.
"""
function graphene_generator(L::Integer, C::Integer, CC_bond::Float64 = 1.422, bc::String = "obc")
    if L < 4
        throw(DomainError("L has to be at least 4"))
    elseif C < 2
        throw(DomainError("C has to be at least 2"))
    elseif L%2 == 1
        throw(DomainError("L has to be even"))
    elseif bc != "obc"
        throw(DomainError("Currently only obc are supported"))
    end

    xproj  = CC_bond * cos(pi/6)            # Projection of carbon atom on x-axis
    yproj  = CC_bond * cos(pi/3)            # Projection of carbon atom on y-axis

    coords      = zeros(L*C, 2)             # Initialize matrix for coords, only 2 dimensions
    coords[:,1] = repmat(collect(0:xproj:xproj*(L-1)), C, 1)        # Putting x coordinates

    ytmp0 = repmat([yproj;0], round(Int, L/2), 1)                   # Temporary y coords. of C atoms
    ytmp1 = repmat([yproj+CC_bond;2*CC_bond], round(Int, L/2), 1)   # Two different y coords for each polyene chain

    # Looping over the number of polyene chains to be stacked on top of each other
    for j = 1:C       # NOTE: should add polyene chains one by one, is it a problem?
        if j%2 == 1
            coords[(j-1)*L+1:j*L,2] = ytmp0+(j-1)*3/2*CC_bond       # One of the two heights of C atoms
        elseif j%2 == 0
            coords[(j-1)*L+1:j*L,2] = ytmp1+(j-2)*3/2*CC_bond       # The second height
        end
    end

    ## Construct connections
    dist_mat = ones(L*C,L*C)*(2*L*C)        # Initialize distance matrix
    bonds = Set{Tuple}()                    # Initialize bonds NOTE: deprecated
    for i = 1:L*C
        dist_mat[i,i] = 0
    end
    for j = 1:C     # Loop over polyene chains
        for i = 1:L     # Loop over atoms in the chain
            if i < L
                push!(bonds,(i+(j-1)*L,i+1+(j-1)*L))
                push!(bonds,(i+1+(j-1)*L,i+(j-1)*L))
                dist_mat[i+(j-1)*L,i+1+(j-1)*L] = 1
                dist_mat[i+1+(j-1)*L,i+(j-1)*L] = 1
            end
            if j%2 == 1 && i%2 == 1 && i < L && j < C
                push!(bonds,(i+(j-1)*L,(i+j*L)%(L*C)))
                push!(bonds,((i+j*L)%(L*C),i+(j-1)*L))
                dist_mat[i+(j-1)*L,(i+j*L)%(L*C)] = 1
                dist_mat[(i+j*L)%(L*C),i+(j-1)*L] = 1
            elseif j%2 == 0 && i%2 == 0 && j < C
                push!(bonds,(i+(j-1)*L,(i+j*L)))
                push!(bonds,((i+j*L),i+(j-1)*L))
                dist_mat[i+(j-1)*L,(i+j*L)] = 1
                dist_mat[(i+j*L),i+(j-1)*L] = 1
            end
        end
    end

    # Compute all the distances between atoms
    populate_dist_mat!(dist_mat)

    return Topology(2, L*C, coords, bonds, dist_mat, bc)
end


"""
    nanotube_generator(n, m , l, CC_bond = 1.422, bc = "obc")

Generate a the topology of a (n,m) carbon nanotube of length l.
"""
function nanotube_generator(n::Integer, m::Integer, l::Integer, CC_bond::Float64 = 1.422, bc::String = "obc")
    if n < m || n+m < 6
        throw(DomainError("n + m must be at least 6, with n >= m"))
    elseif bc != "obc"
        throw(DomainError("Currently only obc are supported"))
    end

    if n == m
        L = (l+1)*2     # Polyene chain length
        C = n+m+2       # Number of polyene chains
    elseif m == 0
        L = n*2
        C = l+1
    else
        throw(ErrorException("Chiral nanotubes are not supported yet!"))
    end

    graphene = graphene_generator(L, C, CC_bond)    # Generate a graphene nanoribbon
    bonds    = graphene.bonds                       # Inherit the bonds
    dist_mat = graphene.dist_mat                    # Inherit the distance matrix
    coords   = zeros(L*C, 3)                        # This is required to avoid a scope error

    if n == m           # Armchair nanotube
        z   = graphene.coords[:,1]              # The x-axis of the graphene ribbon becomes the z-axis
        xgr = graphene.coords[:,2]              # While the y-axis becomes the x-axis which will be wrapped

        # Compute the radius of the nanotube
        sinfac = (maximum(xgr)+CC_bond)/(2*pi)  # Sinusoidal factor
		alpha  = CC_bond/sinfac                 # Angle
		radius = (CC_bond/2)/(sin(alpha/2))     # Radius of the nanotube

        # Roll x- and y-axis of the graphene ribbon
        x = radius * cos(1/sinfac * xgr)
        y = radius * sin(1/sinfac * xgr)

        # Shift the nanotube such that is centered along the z-axis
		zshift = abs(maximum(z)-minimum(z))/2
		z = z - zshift

        # Connect the edges of the ribbon by adding the bonds to the list and modifying the distance matrix
        for i = 1:L
            if i%2 == 0
                push!(bonds,(i,i+(C-1)*L))
                push!(bonds,(i+(C-1)*L,i))
                dist_mat[i,i+(C-1)*L] = 1
                dist_mat[i+(C-1)*L,i] = 1
            end
        end
        populate_dist_mat!(dist_mat)               # Update distance matrix
        coords = [x y z]                           # Bind coordinates in a single matrix
    elseif m == 0       # Zig zag nanotube
        xgr = graphene.coords[:,1]                 # The x-axis of the graphene ribbon which will be wrapped
        z   = graphene.coords[:,2]                 # The y-axis of the graphene ribbon becomes the z-axis

        # Compute the radius of the nanotube
		xproj_CC_bond = abs(xgr[2] - xgr[1])            # Carbon-carbon bond length projected onto the x-axis
        sinfac = (maximum(xgr)+ xproj_CC_bond)/(2*pi)   # Sinusoidal factor
		radius = (xproj_CC_bond/2)/(sin(pi/L))          # Radius of the nanotube

        # Roll x- and y-axis of the graphene ribbon
        x = radius * cos(1/sinfac * xgr)
        y = radius * sin(1/sinfac * xgr)

        # Shift the nanotube such that is centered along the z-axis
		zshift = abs(maximum(z)-minimum(z))/2
		z = z - zshift

        # Connect the edges of the ribbon by adding the bonds to the list and modifying the distance matrix
        for j = 1:C
                push!(bonds,(1+(j-1)*L,j*L))
                push!(bonds,(j*L,1+(j-1)*L))
                dist_mat[1+(j-1)*L,j*L] = 1
                dist_mat[j*L,1+(j-1)*L] = 1
        end
        populate_dist_mat!(dist_mat)             # Update distance matrix
        coords = [x y z]                         # Bind coordinates in a single matrix
    else
        throw(ErrorException("I don't know how you got here!"))
    end

    return Topology(3, L*C, coords, bonds, dist_mat, bc)
end


"""
    populate_dist_mat!(D)

Populate the matrix `D` with the shortest path length between all the vertex of a graph in O(N^3 + logN) complexity.

# Arguments
* `D::Array{Int64,2}`: A matrix containing zeros on the diagonal, ones on connecting vertex and all the rest 2*dim(D,1).
"""
function populate_dist_mat!(D)
    L = size(D,1)
    Dold = deepcopy(D)
    for iter = 1:L
        for i = 1:L
            for j = i+1:L
                D[i,j] = D[j,i] = min(D[i,j],minimum(D[i,:]+D[:,j]))
            end
        end
        if norm(Dold-D) == 0
            return
        end
        Dold = deepcopy(D)
    end
end
