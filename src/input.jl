# Input interface for JuHam

# Wrap all input parameters such that we can run in parallel many inputs using pmap
type Input
    topology::Topology
    hamiltonian::Hamiltonian
    parameters::Dict
end

function polyene_model_inpgen(N::Int64, eta::Real, operator::Function, pbc::Bool)
    topo = one_dim_lattice_generator(N,1,pbc)
    ham = dim_tb_hamiltonian(topo, eta, pbc)
    dic = Dict("eta"=>eta, "N"=>N, "tps_op"=>operator, "pbc"=>pbc)

    ret = Input(topo,ham,dic)
    return ret
end

function polyene_real_inpgen(N::Int64, eta::Real, operator::Function, pbc::Bool)
    topo = polyene_generator(N,1,2,pbc)
    ham = dim_tb_hamiltonian(topo, eta, pbc)
    dic = Dict("eta"=>eta, "N"=>N, "tps_op"=>operator, "pbc"=>pbc)

    ret = Input(topo,ham,dic)
    return ret
end

function polyene_circle_inpgen(N::Int64, eta::Real, operator::Function, pbc::Bool)
    topo = circled_polyene_generator(N,1)
    ham = dim_tb_hamiltonian(topo, eta, true)
    dic = Dict("eta"=>eta, "N"=>N, "tps_op"=>operator, "pbc"=>true)

    ret = Input(topo,ham,dic)
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
