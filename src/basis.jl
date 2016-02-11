abstract Basis

type DeltaDirac <: Basis
    dim::Integer            # Dimension of the basis
    centers::Array          # Array containing the cartesian coordinates of the sites
end
