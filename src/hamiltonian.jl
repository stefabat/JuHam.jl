
"""
    AbstractHamiltonian

Abstract supertype for all Hamiltonians.
"""
abstract type AbstractHamiltonian end

"""
    Hamiltonian <: AbstractHamiltonian

Define a Hamiltonian for a given physical `model` and `topology`.
The Hamiltonian is expressed in the given `basis` and stored in `matrix`.
"""
struct Hamiltonian <: AbstractHamiltonian
    model   ::Model                 # The model used to construct the Hamiltonian matrix
    topology::Topology              # The topology of the system studied
    basis   ::AbstractBasis         # The basis in which the Hamiltonian is expressed
    matrix  ::Matrix{AbstractFloat} # The actual matrix holding the data
end


"""
    Hamiltonian(mod::Huckel, top::Chain)

Hückel Hamiltonian for a 1D chain.
"""
function Hamiltonian(mod::Huckel, top::Chain)
    mat = SymTridiagonal(fill(mod.alpha,top.L),fill(mod.beta,top.L-1))
    return Hamiltonian(mod, top, Basis(), mat)
end


"""
    Hamiltonian(mod::Huckel, top::Molecule)

Hückel Hamiltonian for an arbitrary molecule.

Note that only carbon atoms are considered as interacting sites.
"""
function Hamiltonian(mod::Huckel, top::Molecule)
    n = count(el->isequal(el,"C"),top.types)
    idx = find(el->isequal(el,"C"),top.types)
    mol = Molecule(n, top.types[idx], top.coords[idx,:])
    mat = diagm(fill(mod.alpha,n))
    for key in keys(mol.bonds)
        mat[key[1],key[2]] = mat[key[2],key[1]] = mod.beta
    end
    assert(issymmetric(mat))
    return Hamiltonian(mod, mol, Basis(), mat)
end

# Generate dimerized Hueckel Hamiltonian for a linear polyene chain
# NOTE: this is a very specific generation of the Hamiltonian and should be used carefully
function generate_hamiltonian(model::DimHuckel, topology::Topology, basis::Basis)
    L = basis.size                          # Basis dimension
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

# Wrapper function to compute the orbital energies and the wavefunction
function solve(H::Hamiltonian)
    # Exact diagonalization
    return eig(H.matrix)
end
