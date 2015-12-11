# Input interface for JuHam

# Wrap all input parameters such that we can run in parallel many inputs using pmap
type Input
    topology::Topology
    hamiltonian::Hamiltonian
    parameters::Dict
	meas_op::Function
end

function polyene_model_inpgen(N::Int64, eta::Real, operator::Function, pbc::Bool)
    top = one_dim_lattice_generator(N,1,pbc)
    ham = dim_tb_hamiltonian(top, eta, pbc)
    dic = Dict("eta"=>eta, "N"=>N, "pbc"=>pbc)
	op  = operator

    ret = Input(top,ham,dic,op)
    return ret
end

function polyene_real_inpgen(N::Int64, eta::Real, operator::Function, pbc::Bool)
    top = polyene_generator(N,1,2,pbc)
    ham = dim_tb_hamiltonian(top, eta, pbc)
    dic = Dict("eta"=>eta, "N"=>N, "pbc"=>pbc)
	op  = operator

    ret = Input(top,ham,dic,op)
    return ret
end

function polyene_circle_inpgen(N::Int64, eta::Real, operator::Function, pbc::Bool)
    top = circled_polyene_generator(N,1)
    ham = dim_tb_hamiltonian(top, eta, true)
    dic = Dict("eta"=>eta, "N"=>N, "pbc"=>true)
	op  = operator

    ret = Input(top,ham,dic,op)
    return ret
end

function nanotube_inpgen(n::Int64, m::Int64, l::Int64, operator::Function, pbc::Bool)
	top = nanotube_generator(n,m,l)
	ham = tb_hamiltonian(top)
	dic = Dict("n"=>n, "m"=>m, "l"=>l, "pbc"=>pbc)
	op  = operator

	ret = Input(top,ham,dic,op)
	return ret
end

function generate_multiple_inputs(generator::Function, Nval, etaval, operator, pbc)
    inplist = Array(Input,0)
    if pbc
        for N in Nval
            for eta_ in etaval
                push!(inplist,generator(N,eta_,x->N/(2*pi)*sin(2*pi/N*x),pbc)) # without phase because x starts from -N/2
            end
        end
    else
        for N in Nval
            for eta_ in etaval
                push!(inplist,generator(N,eta_,operator,pbc))
            end
        end
    end
    return inplist
end
