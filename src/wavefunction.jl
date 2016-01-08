# Wave function class
type WaveFunction
	basis   ::Basis                       # Atomic/site basis of the wave function
    MOcoeffs::Array{Float64,2}            # Molecular orbital coefficeints
end

