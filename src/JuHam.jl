
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

include("wavefunction.jl")
export RBM, get_Nh, get_Nσ, log_psi, log_psi_over_psi, psi_over_psi, update_θ


export readxyz
# export DeterminantBasis
export Huckel
export Hamiltonian
export solve

end
