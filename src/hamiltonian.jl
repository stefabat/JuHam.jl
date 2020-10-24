

"""
    AbstractHamiltonian

Abstract supertype for all Hamiltonians.
"""
abstract type AbstractHamiltonian end


"""
    Hamiltonian <: AbstractHamiltonian

Represent a `Model` Hamiltonian for a system with a given `Topology`.
"""
struct Hamiltonian <: AbstractHamiltonian
    model   ::Model                     # The model used to construct the Hamiltonian matrix
    topology::Topology                  # The topology of the system studied
    matrix  ::Matrix{Float64}           # The actual matrix holding the data
end


"""
    Hamiltonian(model::TightBinding, topology::Chain)

Construct the tight-binding Hamiltonian for a one-dimensional `Chain`.
"""
function Hamiltonian(model::TightBinding, topology::Chain)
    matrix = SymTridiagonal(fill(model.α,topology.L),fill(model.β,topology.L-1))
    return Hamiltonian(model, topology, matrix)
end


"""
    Hamiltonian(model::TightBinding, molecule::Molecule)

Construct the tight-binding Hamiltonian for a `Molecule`.

All non-hydrogen atoms are considered as interacting centers.
"""
function Hamiltonian(model::TightBinding, molecule::Molecule)
    # identify the indices of non-hydrogen atoms
    ids = findall(!=("H"),molecule.atoms)
    M = size(ids,1)
    # create an MxM matrix with α on the diagonal
    matrix = diagm(fill(model.α, M))

    for j=1:M
        for i=j+1:M
            if molecule.D[ids[i],ids[j]] == 1
                matrix[i,j] = matrix[j,i] = model.β
            end
        end
    end

    return Hamiltonian(model, molecule, matrix)
end



"""
    solve(H::Hamiltonian)

Compute eigenvalues and eigenvectors of the `Hamiltonian` matrix.
"""
function solve(H::Hamiltonian)
    # exact diagonalization
    return eigen(H.matrix)
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