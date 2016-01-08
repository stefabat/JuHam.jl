# Hamiltonian composite type
type Hamiltonian
    model   ::Model
    topology::Topology
    basis   ::Basis
    matrix  ::Array{Float64,2}
end

# Generate Huckel Hamiltonian
function generate_hamiltonian(model::HuckelModel, topology::Topology)
    N      = size(topology.coords,1)            # Number of sites
    matrix = diagm(repmat([model.alpha],N))     # Coulumb integral on the diagonal
    for bond in topology.bonds
        matrix[bond[1],bond[2]] = model.beta    # Bond integrals according to topoogy
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

