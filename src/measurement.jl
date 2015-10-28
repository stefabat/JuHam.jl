# Define measurements types and functions

# Measurement type
type measurement
    exp_values::Array{Float64}              # Expectation values
end

# Function to compute the TPS in any dimension with an arbitrary operator
function compute_tps(hamiltonian::Hamiltonian, topology::Topology, operator::Function)
    xyz = operator(topology.xyz)
    eig_vectors = hamiltonian.eig_vectors
    N = hamiltonian.nrows
    n_electrons = convert(Int64,N/2)
    tps = zeros(1,3)
    ndim = topology.ndim
    for i = 1:n_electrons
        for j = n_electrons+1:N
            eigensquared = eig_vectors[:,i].*eig_vectors[:,j]
            for k = 1:ndim
                tps[k] += 2.0*dot(eigensquared,xyz[:,k])^2
            end
        end
    end
    return tps
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
