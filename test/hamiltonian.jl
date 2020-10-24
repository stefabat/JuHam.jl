
@testset "hamiltonian" begin

    chain = Chain(10, 0.5)
    H = Hamiltonian(huckel(-1), chain)
    Evals,Evecs = solve(H)

    # analytic energy for open linear chain
    energy(n::Integer) = -2.0*cos(n*pi/11)
    @test Evals ≈ energy.(1:10)

    benzene = Molecule("benzene.xyz")
    H = Hamiltonian(huckel(-1), benzene)
    Evals,Evecs = solve(H)
    
    # analytic energies for benzene
    Eanalytic = [-2.0,-1.0,-1.0,1.0,1.0,2.0]
    @test Evals ≈ Eanalytic
    
end