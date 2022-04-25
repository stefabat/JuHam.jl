
abstract type Model end

"""
    TightBinding <: Model

The `TightBinding` model with on-site integrals equal to `α` and
hopping integrals equal to `β`.

The on-site integral is α = ⟨ϕₘ|Ĥ|ϕₘ⟩
The hopping integral is β = ⟨ϕₘ|Ĥ|ϕₙ⟩
The overlap integral is γ = ⟨ϕₘ|ϕₙ⟩
"""
struct TightBinding <: Model
    α::Float64      # on-site integral
    β::Float64      # hopping integral
    γ::Float64      # overlap integral
end


"""
    huckel( β::Real)

Construct a `TightBinding` model which corresponds to the Hückel Hamiltonian.

The on-site and overlap integrals `α` are assumed to be zero and only the
hopping integral `β` has to be specified.
"""
function huckel(β::Real)
    return TightBinding(0.0, β, 0.0)
end


struct Ising <: Model
    hfield::Float64
end



# """
#     ExtendedHuckel <: Model

# The extended Hückel model.

#     ExtendedHuckel()

# Create an extended Hückel model.
# """
# struct ExtendedHuckel <: Model
# end


# """
#     Hubbard <: Model

# The Hubbard model.

#     Hubbard(t, U, N, d)

# Create a Hubbard model with hopping `t`, on-site 2e repulsion `U` and `N` electrons.
# `d` is the maximal distance of interaction between sites.
# """
# struct Hubbard <: Model
#     t::AbstractFloat    # Hopping integrals (like in TB)
#     U::AbstractFloat    # On-site repulsion
#     N::Integer          # Number of electrons
#     d::Integer          # Max distance of interaction for tij
# end
