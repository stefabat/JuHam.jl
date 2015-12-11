# Define measurements types and functions

# Measurement type
type measurement
    exp_values::Array{Float64}              # Expectation values
end

# Function to compute the TPS in any dimension with an arbitrary operator
function compute_tps(hamiltonian::Hamiltonian, topology::Topology, operator::Function)
    N = hamiltonian.nrows
    M = convert(Int64,N/2)
    ndim = topology.ndim
    xyz = zeros(N,3)
    for i = 1:ndim
        xyz[:,i] = operator(topology.xyz[:,i])
    end
	eig_vectors = hamiltonian.eig_vectors
    tps = zeros(1,3)
    for i = 1:M
        for j = M+1:N
            eig_squared = eig_vectors[:,i].*eig_vectors[:,j]
            for k = 1:ndim
                tps[k] += 2.0*dot(eig_squared,xyz[:,k])^2
            end
        end
    end
    return tps
end

# Function to compute the average position <Psi|X|Psi> of a state |Psi>
function compute_avgpos(hamiltonian::Hamiltonian, topology::Topology, operator::Function)
	N = hamiltonian.nrows
	M = convert(Int64,N/2)
	ndim = topology.ndim
	xyz = zeros(N,3)
	for i = 1:ndim
		xyz[:,i] = operator(topology.xyz[:,i])
	end
	eig_vectors = hamiltonian.eig_vectors
	avgpos = zeros(1,3)
	for i = 1:M
		eig_squared = eig_vectors[:,i].*eig_vectors[:,i]
		for k = 1:ndim
			avgpos[k] += dot(eig_squared,xyz[:,k])
		end
	end
	return avgpos
end


# Function to compute the polarizability in any dimension with an arbitrary operator
# Note: for the moment restricted to parallel directions of the cartesian components
# TODO: compute also xy, xz and yz polarizabilities
function compute_polarizability(hamiltonian::Hamiltonian, topology::Topology, operator::Function)
    xyz = operator(topology.xyz)
    eig_vectors = hamiltonian.eig_vectors
    eig_values = hamiltonian.eig_values
    N = hamiltonian.nrows
    n_electrons = convert(Int64,N/2)
    polarizability = zeros(1,3)
    ndim = topology.ndim
    for i = 1:n_electrons
        for j = n_electrons+1:N
            eigensquared = eig_vectors[:,i].*eig_vectors[:,j]
            for k = 1:ndim
                polarizability[k] += 4.0*(dot(eigensquared,xyz[:,k])^2)/(eig_values[j] - eig_values[i])
            end
        end
    end
    return polarizability
end

# Parallel version: NOT WORKING
# Function to compute the TPS in any dimension with an arbitrary operator
function compute_tps_parallel(hamiltonian::Hamiltonian, topology::Topology, operator::Function)
    xyz = convert(SharedArray,operator(topology.xyz))
    eig_vectors = convert(SharedArray,hamiltonian.eig_vectors)
    N = hamiltonian.nrows
    n_electrons = convert(Int64,N/2)
    tps = zeros(1,3)
    ndim = topology.ndim
    for coord = 1:ndim
        for i = 1:n_electrons
            tps[coord] = @parallel (+) for j = n_electrons+1:N
                2.0*dot((eig_vectors[:,i].*eig_vectors[:,j]),xyz[:,coord])^2
            end
        end
    end
    return tps
end

function test_paral(N)
    res = zeros(1,2)
    tic()
    for j = 1:2
        for k=0.1:0.1:0.5
            res[j] = @sync @parallel (+) for i=1:N
                k*(-1)^(j-1)
            end
        end
    end
    s=toc();
    println("Number of Heads: $res in $s seconds")
end
