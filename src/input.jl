# Input interface for JuHam

# Wrap all input parameters such that we can run in parallel many inputs using pmap
type Input
    topology::Topology
    hamiltonian::Hamiltonian
    parameters::Dict
end

function model_polyene_inpgen(N::Int64, eta::Real, operator::Function, pbc::Bool)
    topo = one_dim_lattice_generator(N,1,pbc)
    ham = dim_tb_hamiltonian(topo, eta, pbc)
    dic = Dict("eta"=>eta, "N"=>N, "tps_op"=>operator, "pbc"=>pbc)

    ret = Input(topo,ham,dic)
    return ret
end

function real_polyene_inpgen(N::Int64, eta::Real, operator::Function, pbc::Bool)
    topo = polyene_generator(N,1,pbc)
    ham = dim_tb_hamiltonian(topo, eta, pbc)
    dic = Dict("eta"=>eta, "N"=>N, "tps_op"=>operator, "pbc"=>pbc)

    ret = Input(topo,ham,dic)
    return ret
end

function circled_polyene_inpgen(N::Int64, eta::Real, operator::Function)
    topo = circled_polyene_generator(N,1)
    ham = dim_tb_hamiltonian(topo, eta, true)
    dic = Dict("eta"=>eta, "N"=>N, "tps_op"=>operator)

    ret = Input(topo,ham,dic)
    return ret
end
