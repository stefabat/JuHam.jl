# Simulation function

function simulation(input::Input)
    topology = input.topology
    hamiltonian = input.hamiltonian
    parameters = input.parameters
    tps_op = input.meas_op
    pbc =  parameters["pbc"]

    diagonalize_hamiltonian!(hamiltonian)
    tps = compute_tps(hamiltonian, topology, tps_op)
	println("<Psi|X^2|Psi>: ",tps[1])
	avgpos = compute_avgpos(hamiltonian, topology, tps_op)
	println("<Psi|X|Psi>: ",avgpos[1])

    ret = Dict("topology"=>topology, "energies"=>hamiltonian.eig_values, "parameters"=>parameters,"tps"=>tps,"exec_time"=>t)
    return ret
end

