# Simple scripting page

for n = 250:500:1750
	for shift = 0:1.1:10
		tic();
		println("###########################################")
		println("N: ",n,"\tshift: ",shift)
		input = polyene_model_inpgen(n,1.1,x->n/(2*pi)*sin(2*pi/n*(x+shift)),true)
		#= input = polyene_model_inpgen(n,1.0,x->x+shift,true) =#
		topology = input.topology
		hamiltonian = input.hamiltonian
		parameters = input.parameters
		tps_op = input.meas_op
		pbc =  parameters["pbc"]

		JuHam.diagonalize_hamiltonian!(hamiltonian)
		tps = JuHam.compute_tps(hamiltonian, topology, tps_op)
		@printf "<Psi|(X+%.2f" shift
		println(")^2|Psi>:\t",tps[1])
		avgpos = JuHam.compute_avgpos(hamiltonian, topology, tps_op)
		@printf "<Psi|X+%.2f" shift
		println("|Psi>^2:\t",avgpos[1]^2)
		t=toc();
		println("###########################################\n")
	end
end

#ret = Dict("topology"=>topology, "energies"=>hamiltonian.eig_values, "parameters"=>parameters,"tps"=>tps,"exec_time"=>t)
