## Topology objects

# Base composite type describing a topology
type Topology
    xyz::Array{Float64}         # N x ndim matrix containing the x, y and z coordinates in the 1st, 2nd and 3rd column, respectively
    ndim::Int64                 # Number of active dimensions
    bonds::Array{Tuple}         # Array of tuples. Each tuple describes a bond between two atoms

    # Simplified constructor if no bonds exist between the atoms
    function Topology(xyz, ndim)
        bonds = Array(Tuple,0)
        new(xyz, bonds)
    end
    # Default constructor with bonds (explicitly needed because of simplified constructor above)
    function Topology(xyz, ndim, bonds)
        new(xyz, ndim, bonds)
    end
end

# Generate the topology of a polyene chain of length nsites
function polyene_generator(nsites, bond_ratio, two_bond_length)
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
    xyz[:,2] = flipdim(xyz[:,2],1)                                 # Flip elements such that first bond goes downward (convention)
    xyz[:,1] = xyz[:,1] - xproj*(nsites-1)/2                    # Shift molecules, such that origin is in the middle of the chain

    ## Construct connections
    bonds = Array(Tuple,0)
    for i = 1:nsites-1                                          # Push bonds into list, no PBC in this case
        push!(bonds,(i,i+1))
        push!(bonds,(i+1,i))
    end

    ## Generate Topology
    ret = Topology(xyz, 2, bonds)
    return ret
end

# Generate the topology of a graphene sheet of with C polyene chaines, each of length N
function graphene_generator(N,C)
    CC_bond = 1                             # Length of C-C bond
    xproj = CC_bond * cos(pi/6)             # Projection of carbon atom on x-axis
    zproj = CC_bond * cos(pi/3)             # Projection of carbon atom on z-axis
    xyz = zeros(N*C,3)
    xyz[:,1] = repmat(collect(0:xproj:xproj*(N-1)),C,1)
    ztmp0 = repmat([zproj;0],convert(Int64,N/2),1)
    ztmp1 = repmat([zproj+CC_bond;2*CC_bond],convert(Int64,N/2),1)
    for j = 1:C       # ERROR: should add polyene chains one by one
        if j%2 == 1
            xyz[(j-1)*N+1:j*N,3] = ztmp0+(j-1)*3/2*CC_bond
        elseif j%2 == 0
            xyz[(j-1)*N+1:j*N,3] = ztmp1+(j-2)*3/2*CC_bond
        end
    end

    ## Construct connections
    bonds = Array(Tuple,0)
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

    ## Generate Topology
    ret = Topology(xyz, 3, bonds)
    return ret
end

# Generate the topology of a (n,m) nanotube of length l
function nanotube_generator(n,m,l)
    assert(n+m>3 || n>4)
    if n == m
        println("Armchair nanotube")
        N = (l+1)*2     # Polyene chain length
        C = n+m+2       # Number of polyene chains
    elseif m == 0
        println("Zigzag nanotube")
        N = n*2
        C = l+1
    else
        throw("Chiral nanotubes not supported yet!")
    end
    graphene = graphene_generator(N,C)
    bonds = graphene.bonds
    #println("graphene bonds: ", size(bonds))
    xyz = graphene.xyz
    if n == m               # Armchair nanotube
        z = graphene.xyz[:,1]           # Note that the x-axis of graphene becomes the z-axis
        xgr = graphene.xyz[:,3]         # And viceversa
        CC_bond = xgr[N+1]-xgr[1]
        fac = (maximum(xgr)+CC_bond)/(2*pi)         # Compute factor for nanotube radius
        x = fac * cos(1/fac * xgr)                  # Roll x- and y-axis of graphene sheet
        y = fac * sin(1/fac * xgr)
        for i = 1:N                     # Add new bonds connecting the edgese to the list
            if i%2 == 0
                push!(bonds,(i,i+(C-1)*N))
                push!(bonds,(i+(C-1)*N,i))
            end
        end
        xyz = [x y z]
    elseif m == 0
        xgr = graphene.xyz[:,1]
        z = graphene.xyz[:,3]
        CC_bond = z[N+1]-z[1]
        fac = (maximum(xgr)+CC_bond)/(2*pi)
        x = fac * cos(1/fac * xgr)
        y = fac * sin(1/fac * xgr)
        for j = 1:C
                push!(bonds,(1+(j-1)*N,j*N))
                push!(bonds,(j*N,1+(j-1)*N))
        end
        xyz = [x y z]
    else
        throw("I don't know you got here!")
    end

    ret = Topology([x y z], 3, bonds)      # Create return Topology object
    return ret
end

# Simple helper function to visualize any topology
function plot_topology(topology)
    x = topology.xyz[:,1]
    y = topology.xyz[:,2]
    z = topology.xyz[:,3]
    edges = topology.bonds
    hold(true)
    for i in edges
        plot3D([x[i[1]],x[i[2]]],[y[i[1]],y[i[2]]],[z[i[1]],z[i[2]]],"-*k")
    end
    axis(xmin=minimum(x)-1,xmax=maximum(x)+1,ymin=minimum(y)-1,ymax=maximum(y)+1,zmin=minimum(z)-1,zmax=maximum(z)+1)
    hold(false)
end
