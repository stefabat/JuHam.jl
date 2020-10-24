
using Test
using JuHam

tests = ["molecule","hamiltonian"]

for test in tests
    include("$test.jl")
end
