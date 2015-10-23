# Define parametric abstract type for Hamiltonians
# T defines the type of the matrix
abstract Hamiltonian{T}

type TightBinding{Real} <: Hamiltonian{Real}
    matrix::Array{T,2}          # Hamiltonian matrix always of dimension 2
    nrwos::Int64                # Int64 according to output of size() function
    ncols::Int64                # Uint64 would probably be better?
    eta::Float64                # beta1/beta2
    sites::Int64                # Number of sites
    nelectrons::Int64           # Number of electrons
    pbc::bool                   # Periodic boundary condition
end
