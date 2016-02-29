# The Model type defines the model Hamiltonian
# Each model contains the parameters characterizing it
abstract Model

# The simple Hueckel model
# Only two parameters, alpha and beta
type HuckelModel <: Model
    alpha::Float64      # Coulomb integral
    beta ::Float64      # Resonance/bond integral
    max_dist::Integer   # Max distance of interaction
end

# The dimerized simple Hueckel model
# Assuming eta = beta1/beta2 and beta1+beta2=const
type DimHuckelModel <: Model
    alpha::Float64      # Coulomb integral
    eta  ::Float64      # beta1/beta2 ratio
    max_dist::Integer   # Max distance of interaction
end

# The tight-binding model
# The tij parameter is a function of the site indeces i and j
type TightBinding <: Model
    tij::Function       # Hopping integrals
    max_dist::Integer   # Max distance of interaction
end
