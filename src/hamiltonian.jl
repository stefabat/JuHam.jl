# Define parametric composite type for Hamiltonians
# T defines the type of the matrix elements
type Hamiltonian
    matrix::Array{Float64,2}          # Hamiltonian matrix always of dimension 2
    nrows::Int64                # Int64 according to output of size() function
    ncols::Int64                # Uint64 would probably be better? Redundant, Hamiltonian always sqaure
    eig_values::Array{Float64}     # Eigenvalues of the Hamiltonian
    eig_vectors::Array{Float64}       # Eigenvectors of the Hamiltonian

    function Hamiltonian(matrix, nrows, ncols)
        eig_values = Array(Float64,0)
        eig_vectors = Array(Float64,0)
        new(matrix, nrows, ncols, eig_values, eig_vectors)
    end
end

# Generate the tight-binding Hamiltonian for a given topology
function tb_hamiltonian(topology)
    xyz = topology.xyz
    bonds = topology.bonds
    N = size(xyz,1)
    matrix = zeros(N,N)
    for bond in bonds
        matrix[bond[1],bond[2]] = 1
    end

    ret = Hamiltonian(matrix, N, N)
    return ret
end

# Wrapper function to compute eigenvalues and eigenvectors of an Hamiltonian object
function diagonalize_hamiltonian(hamiltonian)
    hamiltonian.eig_values,hamiltonian.eig_vectors = eig(hamiltonian.matrix)
    return hamiltonian
end
