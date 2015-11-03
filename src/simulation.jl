# Simulation function

function polyene_tps_simulation(input::Input)
    tic();
    topology = input.topology
    hamiltonian = input.hamiltonian
    parameters = input.parameters
    tps_op = parameters["tps_op"]
    pbc =  parameters["pbc"]

    hamiltonian = diagonalize_hamiltonian(hamiltonian)
    tps = compute_tps(hamiltonian, topology, tps_op)
    polarizability = compute_polarizability(hamiltonian, topology, tps_op)
    t=toc();

    ret = Dict("topology"=>topology, "hamiltonian"=>hamiltonian, "parameters"=>parameters,"tps"=>tps,"polarizability"=>polarizability,"exec_time"=>t)
    return ret
end
