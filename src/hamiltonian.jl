# Parent type of all Hamiltonians
abstract type Hamiltonian

# Tight-Binding Hamiltonian
type TBHamiltonian <: Hamiltonian
    matrix::Array{Float64,2}          # Hamiltonian matrix 
	pbc::Bool




    function Hamiltonian(matrix::Array{Float64,2}, nrows::Int64, ncols::Int64)
        eig_values = Array(Float64,0)
        eig_vectors = Array(Float64,0)
        new(matrix, nrows, ncols, eig_values, eig_vectors)
    end
end

# Generate the tight-binding Hamiltonian for a given topology
function tb_hamiltonian(topology::Topology)
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

# Generate a tight-binding Hamiltonian for two different bond interactions
# NOTE: ad-hoc solution, not general
function dim_tb_hamiltonian(topology::Topology, eta::Float64, pbc::Bool)
    xyz = topology.xyz
    N = size(xyz,1)
    assert(N%2==0)
    matrix = full(SymTridiagonal(fill(0,N),repmat([eta,1/eta],convert(Int64,N/2))))
    if pbc
        matrix[1,end] = matrix[end,1] = 1/eta
    end

    ret = Hamiltonian(matrix, N, N)
    return ret
end

# Wrapper function to compute eigenvalues and eigenvectors of an Hamiltonian object
function diagonalize_hamiltonian!(hamiltonian::Hamiltonian)
    hamiltonian.eig_values,hamiltonian.eig_vectors = eig(hamiltonian.matrix)
    #return hamiltonian
end
