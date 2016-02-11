## Topology objects

# Base composite type describing a topology
type Topology
    dim     ::Integer             # Number of active dimensions (e.g., 1D, 2D, etc.)
    coords  ::Array{Float64}      # N x ndim matrix containing the coordinates of the system
    bonds   ::Set{Tuple}          # Set of tuples. Each tuple describes a bond between two atoms
    dist_mat::Array{Int64,2}      # Distance matrix
end

# Generate the topology of a one dimensional lattice
function one_dim_lattice_generator(L::Integer, lat_const::Real = 1.0)
    coords = collect(linspace(-(L-1)*lat_const/2,(L-1)*lat_const/2,L))  # Generate coordinates of the lattice
    bonds = Set{Tuple}()                # Initialize bonds WARNING: deprecated, soon to be erased
    dist_mat = ones(L,L)*(L+1)          # Initiliaze distance matrix with all entries to L+1
                                        # For a lattice of L sites, L+1 can never be reached
    # Setting nearest neighbors
    for i = 1:L-1
        push!(bonds,(i,i+1))
        push!(bonds,(i+1,i))
        dist_mat[i,i+1] = 1
        dist_mat[i+1,i] = 1
        dist_mat[i,i]   = 0
    end
    dist_mat[L,L] = 0

    # Compute all the distances between the sites
    populate_dist_mat!(dist_mat)

    return Topology(1, coords, bonds, dist_mat)
end

# Generate the topology of a regular polygon
function polygon_generator(N::Integer, bond_length::Real = 1.0)
    ## Construct geometry
    radius = 1 / (2*sin(pi/N))
    coords = zeros(N,2)
    coords[:,1] = radius * cos((2*pi/N)*collect(0:(N-1)))
    coords[:,2] = radius * sin((2*pi/N)*collect(0:(N-1)))

    ## Construct connections
    bonds = Set{Tuple}()
    for i = 1:N-1               # Push bonds into list, no PBC in this case
        push!(bonds,(i,i%N+1))
        push!(bonds,(i%N+1,i))
    end

    return Topology(2, coords, bonds, "obc")
end

# Generate the topology of a graphene ribbon
# NOTE: for the time being, the graphene nanoribbon is
#       generated by stacking polyene chains on top of
#       each other. This might not be the best solution
function graphene_generator(N::Integer, C::Integer, CC_bond::Float64 = 1.0, bc::AbstractString = "obc")
    if N%2 == 1
        throw(ArgumentError("N has to be an even number"))
    elseif C < 2
        throw(ArgumentError("C has to be at least 2"))
    elseif bc != "obc"
        throw(ErrorException("Currently supported boundary conditions are obc"))
    end

    xproj  = CC_bond * cos(pi/6)             # Projection of carbon atom on x-axis
    yproj  = CC_bond * cos(pi/3)             # Projection of carbon atom on y-axis

    coords      = zeros(N*C, 2)
    coords[:,1] = repmat(collect(0:xproj:xproj*(N-1)), C, 1)

    ytmp0 = repmat([yproj;0], round(Int, N/2), 1)
    ytmp1 = repmat([yproj+CC_bond;2*CC_bond], round(Int, N/2), 1)

    for j = 1:C       # ERROR: should add polyene chains one by one
        if j%2 == 1
            coords[(j-1)*N+1:j*N,2] = ytmp0+(j-1)*3/2*CC_bond
        elseif j%2 == 0
            coords[(j-1)*N+1:j*N,2] = ytmp1+(j-2)*3/2*CC_bond
        end
    end

    ## Construct connections
    bonds = Set{Tuple}()
    for j = 1:C
        for i = 1:N
            if i < N
                push!(bonds,(i+(j-1)*N,i+1+(j-1)*N))
                push!(bonds,(i+1+(j-1)*N,i+(j-1)*N))
            end
            if j%2 == 1 && i%2 == 1 && i < N && j < C
                push!(bonds,(i+(j-1)*N,(i+j*N)%(N*C)))
                push!(bonds,((i+j*N)%(N*C),i+(j-1)*N))
            elseif j%2 == 0 && i%2 == 0 && j < C
                push!(bonds,(i+(j-1)*N,(i+j*N)))
                push!(bonds,((i+j*N),i+(j-1)*N))
            end
        end
    end

    return Topology(2, coords, bonds, bc)
end

# Generate the topology of an (n,m) carbon nanotube of length l
function nanotube_generator(n::Integer, m::Integer , l::Integer, CC_bond::Float64 = 1.0, bc::AbstractString = "obc")
    if n < m || n+m < 6
        throw(ArgumentError("n + m must be at least 6, with n >= m"))
    end

    if n == m
        N = (l+1)*2     # Polyene chain length
        C = n+m+2       # Number of polyene chains
    elseif m == 0
        N = n*2
        C = l+1
    else
        throw(ErrorException("Chiral nanotubes are not supported yet!"))
    end

    graphene = graphene_generator(N, C, CC_bond)
    bonds    = graphene.bonds
    coords   = zeros(N*C, 3)

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

        # Connect the edges of the ribbon by adding the bonds to the list
        for i = 1:N
            if i%2 == 0
                push!(bonds,(i,i+(C-1)*N))
                push!(bonds,(i+(C-1)*N,i))
            end
        end
        coords = [x y z]
    elseif m == 0       # Zig zag nanotube
        xgr = graphene.coords[:,1]                 # The x-axis of the graphene ribbon which will be wrapped
        z   = graphene.coords[:,2]                 # The y-axis of the graphene ribbon becomes the z-axis

        # Compute the radius of the nanotube
		xproj_CC_bond = abs(xgr[2] - xgr[1])            # Carbon-carbon bond length projected onto the x-axis
        sinfac = (maximum(xgr)+ xproj_CC_bond)/(2*pi)   # Sinusoidal factor
		radius = (xproj_CC_bond/2)/(sin(pi/N))          # Radius of the nanotube

        # Roll x- and y-axis of the graphene ribbon
        x = radius * cos(1/sinfac * xgr)
        y = radius * sin(1/sinfac * xgr)

        # Shift the nanotube such that is centered along the z-axis
		zshift = abs(maximum(z)-minimum(z))/2
		z = z - zshift

        # Connect the edges of the ribbon by adding the bonds to the list
        for j = 1:C
                push!(bonds,(1+(j-1)*N,j*N))
                push!(bonds,(j*N,1+(j-1)*N))
        end
        coords = [x y z]
    else
        throw(ErrorException("I don't know how you got here!"))
    end

    return Topology(3, coords, bonds, bc)
end

# Populate the distance matrix
# Input: a matrix with zeros on the diagonal and ones on connected
# sites. All the other entries are much larger.
# Ouput: the populated matrix with the distance between each site
# WARNING: very naive algorithm, O(N^3)
function populate_dist_mat!(D)
    L = size(D,1)
    for iter = 1:L-1
        for i = 1:L
            for j = i+2:L
                D[i,j] = minimum(D[i,i+1:L]'+D[i+1:L,j])
                D[j,i] = minimum(D[i,i+1:L]'+D[i+1:L,j])
            end
        end
    end
end

# Generate the topology of a lattice of arbitrary dimension
#=
function lattice_generator(nsites::Array{Integer,1}, lat_consts::Array{Real,1}, bc::String = "obc")
    ndim = length(nsites)       # Number of dimensions
    if ndim != length(lat_consts)
        throw(DimensionMismatch("the arguments nsites and lat_consts have to be one-dimensional arrays of the same length"))
    end
    coords = zeros(nsites,ndim)
    for i = 1:ndim
        coords[:,i] = collect(linspace(-(nsites[i]-1)*lat_consts[i]/2,(nsites[i]-1)*lat_consts[i]/2,nsites[i]))
    end
    bonds = Set{Tuple}()
    for i = 1:ndim
        for j = 1:nsites[i]-1
            push!(bonds,(j + (i-1)*nsites[i]   ,j+1))
        push!(bonds,(i+1,i))
    end
    if pbc
        push!(bonds,(1,N))
        push!(bonds,(N,1))
    end

    # Generate Topology
    ret = Topology(xyz, 1, bonds)
    return ret
end

# Generate the topology of a polyene chain of length nsites
function polyene_generator(nsites, bond_ratio, two_bond_length, pbc)
    ## Construct geometry
    xyz = zeros(nsites,3)
    long_bond = bond_ratio * two_bond_length / (1 + bond_ratio)         # Find the length of the longer bond
    short_bond = two_bond_length - long_bond                            # Find the length of the shroter bond
    xproj = two_bond_length / 2 * cos(pi/6)                             # Projection onto the x-axis of the the second atom
    xtmp = zeros(nsites+1)                                              # Temporary variable
    xtmp[1:2:nsites+1] = 0:2*xproj:(nsites+1)*xproj                     # Fill in with all the x's lying on y=0, one element more
                                                                        # is needed to place correctly the last element
    xyz[1:2:nsites,1] = xtmp[1:2:nsites]                                # Copy over to output, leaving out last element
    for i = 2:2:nsites                  # Find the coordinates for x and y, solving the system of eq. of two circles, bond angle is free
        x = -1/2 * (long_bond^2 - short_bond^2 + xtmp[i-1]^2 - xtmp[i+1]^2) / (xtmp[i+1] - xtmp[i-1])
        y = sqrt(short_bond^2 - x^2 - xtmp[i-1]^2 + 2*x*xtmp[i-1])
        xyz[i,1] = x
        xyz[i,2] = y
    end
    xyz[:,2] = flipdim(xyz[:,2],1) - abs(maximum(xyz[:,2])-minimum(xyz[:,2]))/2                                 # Flip elements such that first bond goes downward (convention)
    xyz[:,1] = xyz[:,1] - xproj*(nsites-1)/2                    # Shift molecules, such that origin is in the middle of the chain

    ## Construct connections
    bonds = Array(Tuple,0)
    for i = 1:nsites-1                                          # Push bonds into list, no PBC in this case
        push!(bonds,(i,i+1))
        push!(bonds,(i+1,i))
    end
    if pbc
        push!(bonds,(1,nsites))
        push!(bonds,(nsites,1))
    end

    ## Generate Topology
    ret = Topology(xyz, 2, bonds)
    return ret
end
=#
