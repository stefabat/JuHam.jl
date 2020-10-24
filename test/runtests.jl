
using Test
using JuHam

tests = ["molecule"]

for test in tests
    include("$test.jl")
end
