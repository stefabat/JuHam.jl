
"""
    AbstractBasis

Abstract supertype for all basis sets.
"""
abstract type AbstractBasis end

struct DummyBasis <: AbstractBasis end

# struct DeltaBasis <: AbstractBasis
#     size::Int
#     ldim::Int
# end

# struct DeterminantBasis <: AbstractBasis
#     alpha::BitMatrix
#     beta ::BitMatrix
# end

# "`N` electrons in `M` MOs"
# function DeterminantBasis(N::Int, M::Int, Ms::Int=1)
#     Nup = Int(N/2)
#     Ndw = Int(N/2)
#     alpha = BitMatrix(binomial(M,Nup),M)
#     beta  = BitMatrix(binomial(M,Ndw),M)
#     i=1
#     for state = 2^M-1:-1:0
#         tmp = count_ones(state)
#         if tmp == Nup && tmp == Ndw
#             str = bin(state,M)
#             alpha[i,:] = BitVector([parse(Bool,str[i]) for i=1:M])
#             beta[i,:] = BitVector([parse(Bool,str[i]) for i=1:M])
#             i += 1
#         elseif tmp == Nup
#             str = bin(state,M)
#             alpha[i,:] = BitVector([parse(Bool,str[i]) for i=1:M])
#             i += 1
#         elseif tmp == Ndw
#             str = bin(state,M)
#             beta[i,:] = BitVector([parse(Bool,str[i]) for i=1:M])
#             i += 1
#         end
#     end
#     return DeterminantBasis(alpha,beta)
# end

# function getsize(basis::DeterminantBasis)
#     return length(basis.alpha)*length(basis.beta)
# end
