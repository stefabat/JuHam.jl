# Wave function class
type WaveFunction
	basis::String							# Name of the basis
    MOenergies::Array{Float64,1}			# Eigenvalues of the Hamiltonian
    MOcoeffs::Array{Float64,2}				# Atomic orbital coefficeintsEigenvectors of the Hamiltonian
end

function basis_dim(wf::WaveFunction)
	return size(wf.MOcoeffs,1)
end

