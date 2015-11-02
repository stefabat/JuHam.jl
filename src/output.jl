## Require HDF5 -- for portability

function format_output(results_arr)
    ntot = size(results_arr,1)
    ret = Array(Dict,ntot)
    for i = 1:ntot
        results = results_arr[i]
        hamiltonian = results["hamiltonian"]
        topology = results["topology"]
        parameters = results["parameters"]
        tps = results["tps"]
        exec_time = results["exec_time"]
        N = parameters["N"]
        eta = parameters["eta"]
        pbc = parameters["pbc"]
        ndim = topology.ndim
        energies = hamiltonian.eig_values
        # orbitals = hamiltonian.eig_vectors
        geometry = topology.xyz[:,1:ndim]
        # tmp = Dict("N"=>N,"eta"=>eta,"pbc"=>pbc,"geometry"=>geometry,"exec_time"=>exec_time,"energies"=>energies,"orbitals"=>orbitals,"tps"=>tps)
        tmp = Dict("N"=>N,"eta"=>eta,"pbc"=>pbc,"geometry"=>geometry,"exec_time"=>exec_time,"energies"=>energies,"tps"=>tps)
        ret[i] = tmp
    end
    return ret
end
