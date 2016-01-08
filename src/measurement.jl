# Measurements function

# Total position spread tensor and polarizability
# of a given wavefunction and topology
# NOTE: only valid in the Huckel approximation
# NOTE: polarizability currently restricted to xx, yy and zz components
# TODO: compute also xy, xz and yz polarizabilities
function compute_tps(wave_function::WaveFunction, topology::Topology, MOenergies)
    N      = wave_function.basis.dim
    M      = round(Int64,N/2)
    ndim   = topology.dim
    coords = zeros(N,ndim)
	eigvec = wave_function.MOcoeffs
    tps    = zeros(1,ndim)
    pola   = zeros(1,ndim)

    # TODO: to parallelize
    for i = 1:M         # Loop over occupied orbitals
        for j = M+1:N   # Loop over empty orbitals
            eigsq = eigvec[:,i].*eigvec[:,j]
            for k = 1:ndim
                dotprod  = dot(eigsq,coords[:,k])^2
                tps[k]  += 2.0*dotprod
                pola[k] += 4.0*dotprod/(MOenergies[j] - MOenergies[i])
            end
        end
    end

    return tps,pola
end

# Function to compute the average position <Psi|X|Psi> of a state |Psi>
function compute_avgpos(wave_function::WaveFunction, topology::Topology)
    N      = wave_function.basis.dim
    M      = round(Int64,N/2)
    ndim   = topology.dim
    coords = zeros(N,ndim)
	eigvec = wave_function.MOcoeffs
	avgpos = zeros(1,ndim)

	for i = 1:M
		eigsq = eigvec[:,i].*eigvec[:,i]
		for k = 1:ndim
			avgpos[k] += dot(eigsq,coords[:,k])
		end
	end

	return avgpos
end

function test_paral(N)
    res = zeros(1,2)
    for j = 1:2
        for k=0.1:0.1:0.5
            res[j] = @sync @parallel (+) for i=1:N
                k*(-1)^(j-1)
            end
        end
    end
end
