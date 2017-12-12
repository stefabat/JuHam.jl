using DataStructures.OrderedDict

"""
    Molecule <: Topology

Composite type describing a molecule.

"""
struct Molecule <: Topology
    natoms::Int
    types ::Array{String}
    coords::Array{AbstractFloat}
    bonds ::OrderedDict{Tuple{Int,Int},AbstractFloat}

    # Constructor
    function Molecule(natoms, types, coords)
    bonds = OrderedDict{Tuple{Int,Int},AbstractFloat}()
    for i = 1:natoms-1
        for j = i+1:natoms
            d = norm(coords[i,:] - coords[j,:])
            if d < cov_radii[types[i]] + cov_radii[types[j]] + 0.25
                push!(bonds,(i,j)=>d)
            end
        end
    end
    new(natoms, types, coords, bonds)
    end

end

"Read an xyz file and generate an object of type `Molecule`"
function readxyz(input)
    f = open(input,"r")
    natoms = parse(readline(f))
    atoms  = readdlm(f, skipstart=1)
    close(f)
    if natoms != size(atoms,1)
        throw(ArgumentError("invalid xyz file: total number of atoms does not match
        the number declared in the first line."))
    end
    types  = convert(Array{String,1}, atoms[:,1])
    coords = convert(Array{AbstractFloat,2}, atoms[:,2:end])
    return Molecule(natoms, types, coords)
end

"Atomic covalent radii taken from http://periodictable.com"
cov_radii = Dict("H" =>0.31,"He"=>0.28,"Li"=>1.28,"Be"=>0.96,"B" =>0.85,
                 "C" =>0.76,"N" =>0.71,"O" =>0.66,"F" =>0.57,"Ne"=>0.58)
