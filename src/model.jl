# Parent type of all models
abstract Model

type HuckelModel <: Model
    alpha::Float64      # Coulomb integral
    beta ::Float64      # Resonance/bond integral
end

