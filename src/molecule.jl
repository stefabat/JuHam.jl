
"""
    Molecule <: Topology

Type representing a molecule consisting of `natoms` atoms of type `types`
"""
struct Molecule <: Topology
    atoms::Vector{String}
    coordinates::Matrix{Float64}
    D::Matrix{Int}      # connectivity matrix
end


"""
    Molecule(atoms::Vector{String}, coordinates::Matrix{Float64}, bond_threshold::Real = 0.25)

Construct a `Molecule` composed by a list of `atoms` and corresponding `coordinates`.

`atoms` is a `Vector{String}` of length N containing the atomic symbols composing the molecule.
Supported atoms are: H, Li, B, C, N, O, F, Na, Al, Si, P, S, Cl, K, Ga, Ge, As, Se and Br.
`coordinates` is a `Matrix{Float64}` of dimension Nx3 containing the coordinates of the atoms.
The coordinates of atom `atoms[I]` are stored in `coordinates[I,:]`.
Two atoms are considered bonded if their distance is less or equal to the sum of their covalent
radii plus `bond_threshold`.
"""
function Molecule(atoms::Vector{String}, coordinates::Matrix{Float64}, bond_threshold::Real = 0.25)
    if size(atoms,1) != size(coordinates,1)
        error("the length of `atoms` and `coordinates` should be equal.")
    end
    N = size(atoms,1)
    D = fill(N,(N,N)) 
    D[1:N+1:end] .= 0
    # determine the connectivity
    for i = 1:N-1
        for j = i+1:N
            d = norm(coordinates[i,:] - coordinates[j,:])
            if d < covalent_radius(atoms[i]) + covalent_radius(atoms[j]) + bond_threshold
                D[i,j] = 1
                D[j,i] = 1
            end
        end
    end
    floydwarshall!(D)

    return Molecule(atoms, coordinates, D)
end


"""
    Molecule(xyzfile::String, bond_threshold::Real = 0.25)

Construct a `Molecule` by reading an xyz file at `xyzfile`.
"""
function Molecule(xyzfile::String, bond_threshold::Real = 0.25)
    atoms,coordinates = readxyz(xyzfile)
    return Molecule(atoms, coordinates, bond_threshold)
end


"""
    readxyz(input::String)

Read an xyz file and return a `Vector{String}` containing the atomic symbols and
a `Matrix{Float64}` containing the atomic coordinates.
"""
function readxyz(input::String)
    f = open(input, "r")
    N = parse(Int, readline(f))
    data = readdlm(f, skipstart=1)
    close(f)
    if N != size(data,1)
        error("invalid xyz file, number of atoms mismatch.")
    end
    atoms = convert(Vector{String}, data[:,1])
    coordinates = convert(Matrix{Float64}, data[:,2:end])

    return atoms,coordinates
end


"""
    covalent_radius(atom::String)

Return the covalent radius of `atom` in Ångstrom.

Data taken from "Cordero et al. Dalt. Trans. 2008, No. 21, 2832.".
"""
function covalent_radius(atom::String)
    atom = uppercase(strip(atom))
    if atom == "H"
        return 0.31
    elseif atom == "LI"
        return 1.28
    elseif atom == "B"
        return 0.84
    elseif atom == "C"
        return 0.76
    elseif atom == "N"
        return 0.71
    elseif atom == "O"
        return 0.66
    elseif atom == "F"
        return 0.57
    elseif atom == "NA"
        return 1.66
    elseif atom == "AL"
        return 1.21
    elseif atom == "SI"
        return 1.11
    elseif atom == "P"
        return 1.07
    elseif atom == "S"
        return 1.05
    elseif atom == "CL"
        return 1.02
    elseif atom == "K"
        return 2.03
    elseif atom == "GA"
        return 1.22
    elseif atom == "GE"
        return 1.20
    elseif atom == "AS"
        return 1.19
    elseif atom == "SE"
        return 1.20
    elseif atom == "BR"
        return 1.20
    else
        error("atom type not supported.")
    end
end
