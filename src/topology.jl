## Topology

type Topology
    bonds::Array{Tuple}
    xyz::Array{Float64}
end

function PolyeneChainGenerator(nsites, bond_ratio, two_bond_length)
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
    xyz[:,2] = flipud(xyz[:,2])                                 # Flip elements such that first bond goes downward (convention)
    xyz[:,1] = xyz[:,1] - xproj*(nsites-1)/2                    # Shift molecules, such that origin is in the middle of the chain

    ## Construct connections
    bonds = Array(Tuple,0)
    for i = 1:nsites-1                                          # Push bonds into list, no PBC in this case
        push!(bonds,(i,i+1))
        push!(bonds,(i+1,i))
    end

    ## Generate Topology
    ret = Topology(bonds, xyz)

    return ret
end

function nanotube_geometry(n,m,l)
    H,N,C = nanotube(n,m,l)                 # Compute tight-binding Hamiltonian
    CC_bond = 1                             # Length of C-C bond
    theta = 2*pi/N                          # Angle between two neighboring carbon atoms
    radius = 1 / (2*sin(theta/2))           # Radius of the nanotube (Assuming C-C bond length is 1)
    vertex_dist = 2                         # Distance of two opposite vertices (in a ring)
    edge_dist = sqrt(3)                     # Distance of two opposite edges (in a ring)
    half_length = vertex_dist*floor(l/2) + vertex_dist/2    # Half length of the entire nanotube (origin z=0)
    xproj = CC_bond * cos(pi/3)         # Projection of carbon atom on x-axis
    xvec = zeros(N*C,1)                     ###
    x0_cell = [xproj,xproj+CC_bond]
    x1_cell = [0,2*CC_bond]
    x0 = repmat(x0_cell,convert(Int64,C/2),1)
    x1 = repmat(x1_cell,convert(Int64,C/2),1)
    for i = 1:C/2
        x0[2*i-1] += (i-1)*3*CC_bond
        x0[2*i] += (i-1)*3*CC_bond
        x1[2*i-1] += (i-1)*3*CC_bond
        x1[2*i] += (i-1)*3*CC_bond
    end
    for i = 1:N/2
        for j = 1:C
            xvec[2*i-1 + (j-1)*N] = x0[j]
            xvec[2*i + (j-1)*N] = x1[j]
        end
    end
    fac = (maximum(xvec)+CC_bond)/(2*pi)
    xcirc = fac * cos(1/fac * xvec)
    ycirc = fac * sin(1/fac * xvec)
    return xcirc,ycirc
end

function graphene(n,m,l)
    H,N,C = nanotube(n,m,l)                 # Compute tight-binding Hamiltonian
    CC_bond = 1                             # Length of C-C bond
    xproj = CC_bond * cos(pi/6)             # Projection of carbon atom on x-axis
    zproj = CC_bond * cos(pi/3)
    xyz = zeros(N*C,3)
    xyz[:,1] = repmat([0:xproj:xproj*(N-1)],C,1)
    ztmp = [repmat([zproj,0],convert(Int64,N/2),1),repmat([zproj+CC_bond,2*CC_bond],convert(Int64,N/2),1)]
    for j = 1:2:C
        xyz[(j-1)*N+1:(j+1)*N,3] = ztmp+(j-1)*3/2*CC_bond
    end

    ## Construct connections
    bonds = Array(Tuple,0)
    for j = 1:C
        for i = 1:N-1
            push!(bonds,(i+(j-1)*N,i+1+(j-1)*N))
            push!(bonds,(i+1+(j-1)*N,i+(j-1)*N))
            if j%2 == 1 && i%2 == 1
                push!(bonds,(i+(j-1)*N,i+j*N))
                push!(bonds,(i+j*N,i+(j-1)*N))
            elseif j%2 == 0 && i%2 == 0 && j < C
                push!(bonds,(i+(j-1)*N,i+j*N))
                push!(bonds,(i+j*N,i+(j-1)*N))
            end
        end
    end

    ## Generate Topology
    ret = Topology(bonds, xyz)

    return ret
end

function plot_geometry(topology)
    x = topology.xyz[:,1]
    z = topology.xyz[:,3]
    edges = topology.bonds
    # Winston.hold(true)
    for i in edges
        Winston.plot([x[i[1]],x[i[2]]],[z[i[1]],z[i[2]]])
    end
    Winston.plot(x,z)
    # hold(false);
end
