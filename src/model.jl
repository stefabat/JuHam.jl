"""
    Model

Abstract supertype for all physical models.
"""
abstract type Model end


"""
    Huckel <: Model

The Hückel model.

    Huckel(alpha::AbstractFloat, beta::AbstractFloat)

Create a Hückel model with the on-site (Coulomb) integral equal to `alpha`
and the hopping (resonance/bond) integral equal to `beta`.
"""
struct Huckel <: Model
    alpha::AbstractFloat      # Coulomb integral
    beta ::AbstractFloat      # Resonance/bond integral
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


# struct PPP <: Model
#     t::AbstractFloat
# end
