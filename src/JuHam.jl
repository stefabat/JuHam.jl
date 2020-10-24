
module JuHam

using LinearAlgebra
using DelimitedFiles
using SparseArrays


include("lattice.jl")
# structs
export Chain

include("molecule.jl")
# structs
export Molecule

include("model.jl")
# structs
export TightBinding
# functions
export huckel

include("hamiltonian.jl")




export readxyz
# export DeterminantBasis
export Huckel
export Hamiltonian
export solve

end
