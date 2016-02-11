# Hamiltonian matrix composite type
# The Hamiltonian matrix for a given model,
# expressed in the basis {\Phi_i}.
# Hij = <\Phi_i| H_eff |\Phi_j>
type Hamiltonian
    model::Model                # The model used to construct the Hamiltonian matrix
    topology::Topology          # The topology of the system studied
    basis::Basis                # The basis in which the Hamiltonian is expressed
    matrix::Array{Float64,2}    # The actual matrix holding the data
    bc::AbstractString          # Boundary conditions
end

# For each model, we have a specialized constructor, which knows the parameters
# of the model chosen.
# Same function name, different signatures: Multiple Dispatch

# Generate Huckel Hamiltonian
function generate_hamiltonian(model::HuckelModel, topology::Topology, basis::Basis)
    L = basis.dim                               # Basis dimension
    matrix = diagm(repmat([model.alpha],N))     # Coulumb integrals on the diagonal
    for bond in topology.bonds
        matrix[bond[1],bond[2]] = model.beta    # Bond integrals according to topology
    end

    return Hamiltonian(model, topology, Basis(N), matrix)
end

# Wrapper function to compute the orbital energies and the wavefunction
function diagonalize_hamiltonian(hamiltonian::Hamiltonian)
    MOcoeffs,MOenergiesc = eig(hamiltonian.matrix)

    return WaveFunction(hamiltonian.basis, MOcoeffs),MOenergies
end

#=
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
=#
