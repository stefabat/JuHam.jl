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


"""
    ExtendedHuckel <: Model

The extended Hückel model.

    ExtendedHuckel()

Create an extended Hückel model.
"""
struct ExtendedHuckel <: Model
end


"""
    Hubbard <: Model

The Hubbard model.

    Hubbard(t, U, n)

Create a Hubbard model with hopping `t` and on-site 2e repulsion `U`.
Sites interacts at a maximal distance `n`.
"""
struct Hubbard <: Model
    t::AbstractFloat    # Hopping integrals (like in TB)
    U::AbstractFloat    # On-site repulsion
    n::Integer          # Max distance of interaction for tij
end


struct PPP <: Model
    t::AbstractFloat
end
# The dimerized simple Hueckel model
# Assuming eta = beta1/beta2 and beta1+beta2=const
struct DimHuckel <: Model
    alpha::Float64      # Coulomb integral
    eta  ::Float64      # beta1/beta2 ratio
end
