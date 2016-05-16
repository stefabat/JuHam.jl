# Hamiltonian matrix composite type
# The Hamiltonian matrix for a given model,
# expressed in the basis {\Phi_i}.
# Hij = <\Phi_i| H_eff |\Phi_j>
type Hamiltonian
    model::Model                # The model used to construct the Hamiltonian matrix
    topology::Topology          # The topology of the system studied
    basis::Basis                # The basis in which the Hamiltonian is expressed
    matrix::Array{Float64,2}    # The actual matrix holding the data
end

# For each model, we have a specialized constructor, which knows the parameters
# of the model chosen.
# Same function name, different signatures: Multiple Dispatch

# Generate simple Hueckel Hamiltonian
function generate_hamiltonian(model::HuckelModel, topology::Topology, basis::Basis)
    L = basis.dim                           # Basis dimension

    # Coulumb integrals on the diagonal, bond integrals on the off diagonals
    matrix = SymTridiagonal(fill(model.alpha,L), fill(model.beta,L-1))

    # If pbc, transform the matrix into a dense one and add the connceting values
    if topology.bc == "pbc"
        matrix = full(matrix)
        matrix[1,end] = matrix[end,1] = model.beta
    end

    return Hamiltonian(model, topology, basis, matrix)
end

# Generate dimerized Hueckel Hamiltonian for a linear polyene chain
# NOTE: this is a very specific generation of the Hamiltonian and should be used carefully
function generate_hamiltonian(model::DimHuckelModel, topology::Topology, basis::Basis)
    L = basis.dim                           # Basis dimension
    # TODO: change it to throw exception
    assert(L%2==0)                          # Check that there are an even number of sites

    # Coulumb integrals on the diagonal and eta, 1/eta on the lower and upper diagonals
    matrix = SymTridiagonal(fill(model.alpha,L), repmat([model.eta,1/model.eta], convert(Int64,L/2)))

    # If pbc, transform the matrix into a dense one and add the connceting values
    if topology.bc == "pbc"
        matrix = full(matrix)
        matrix[1,end] = matrix[end,1] = 1/eta
    end

    return Hamiltonian(model, topology, basis, matrix)
end

# Generate tight-binding Hamiltonian
function generate_hamiltonian(model::TightBinding, topology::Topology, basis::Basis)
    throw(ErrorException("The tight-binding Hamiltonian is not yet implemented"))
    # L = basis.dim                   # Basis dimension
    # matrix = zeros(L, L)            # Initialize data matrix
end

# Wrapper function to compute the orbital energies and the wavefunction
function diagonalize_hamiltonian(hamiltonian::Hamiltonian)
    # Exact diagonalization
    MOcoeffs,MOenergiesc = eig(hamiltonian.matrix)

    return WaveFunction(hamiltonian.basis, MOcoeffs),MOenergies
end
