# Define parametric abstract type for Hamiltonians
# T defines the type of the matrix
type Hamiltonian{T}
    matrix::Array{T,2}        # Hamiltonian matrix always of dimension 2
    nrwos::Int64              # Int64 according to output of size() function
    ncols::Int64              # Uint64 would probably be better?
end

type TightBinding{Real} <: Hamiltonian{Real}
    eta::Float64
    sites::Int64
    pbc::bool
end
