
using LinearAlgebra
using SparseArrays

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
    matrix  ::Matrix{AbstractFloat}       # The actual matrix holding the data
end


"""
    Hamiltonian(mod::Huckel, top::Chain)

Hückel Hamiltonian for a 1D chain.
"""
function Hamiltonian(mod::Huckel, top::Chain)
    mat = SymTridiagonal(fill(mod.alpha,top.L),fill(mod.beta,top.L-1))
    return Hamiltonian(mod, top, DummyBasis(), mat)
end


"""
    Hamiltonian(mod::Huckel, top::Molecule)

Hückel Hamiltonian for an arbitrary molecule.

Note that only carbon atoms are considered as interacting sites.
"""
function Hamiltonian(mod::Huckel, mol::Molecule)
    n = count(el->isequal(el,"C"),mol.types)
    # idx = find(el->isequal(el,"C"),mol.types)
    idx = findall(el->isequal(el,"C"),mol.types)
    top = Molecule(n, mol.types[idx], mol.coords[idx,:])
    mat = diagm(fill(mod.alpha,n))
    for key in keys(top.bonds)
        mat[key[1],key[2]] = mat[key[2],key[1]] = mod.beta
    end
    @assert(issymmetric(mat))
    return Hamiltonian(mod, top, DummyBasis(), mat)
end


# """
#     Hamiltonian(mod::Hubbard, top::Chain)

# Hubbard Hamiltonian for a 1D chain.
# """
# function Hamiltonian(mod::Hubbard, top::Chain)
#     @time basis = DeterminantBasis(mod.N,top.L)
#     Mup = size(basis.alpha,1)
#     Mdw = size(basis.beta,1)
#     println("FCI space dimension: ",getsize(basis))
#     H = spzeros(Mup*Mdw,Mup*Mdw)
#     for bra_dw = 1:Mdw
#         for ket_dw = 1:Mdw
#             d_dw = distance(basis.beta[bra_dw,:],basis.beta[ket_dw,:])
#             if d_dw <= mod.d
#                 for bra_up = 1:Mup
#                     for ket_up = 1:Mup
#                         d_up = distance(basis.alpha[bra_up,:],basis.alpha[ket_up,:])
#                         if d_up <= mod.d
#                             if d_up == d_dw == 0
#                                 H[Mdw*(bra_dw-1)+bra_up,Mdw*(ket_dw-1)+ket_up] = +mod.U
#                             else
#                                 H[Mdw*(bra_dw-1)+bra_up,Mdw*(ket_dw-1)+ket_up] = -mod.t
#                             end
#                         end
#                     end
#                 end
#             end
#         end
#     end
#     return Hamiltonian(mod, top, basis, H)
# end

# Wrapper function to compute the orbital energies and the wavefunction
function solve(H::Hamiltonian)
    # Exact diagonalization
    return eigen(H.matrix)
end

function distance(phi1, phi2)
    # idx = find(xor.(phi1,phi2))
    idx = findall(xor.(phi1,phi2))
    if length(idx) > 0
        return abs(idx[1]-idx[2])
    else
        return 0
    end
end
