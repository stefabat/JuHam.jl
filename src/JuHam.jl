
module JuHam

using LinearAlgebra

include("lattice.jl")
export Chain

include("molecule.jl")
export Molecule

include("model.jl")
include("basis.jl")
# include("wavefunction.jl")
include("hamiltonian.jl")




export readxyz
# export DeterminantBasis
export Huckel
export Hamiltonian
export solve

end
