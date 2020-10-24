
@testset "molecule" begin
  
    benzene = Molecule("benzene.xyz")

    # connectivity matrix of "benzene.xyz"
    D = [0  1  2  3  2  1  1  2  3  4  3  2
         1  0  1  2  3  2  2  1  2  3  4  3
         2  1  0  1  2  3  3  2  1  2  3  4
         3  2  1  0  1  2  4  3  2  1  2  3
         2  3  2  1  0  1  3  4  3  2  1  2
         1  2  3  2  1  0  2  3  4  3  2  1
         1  2  3  4  3  2  0  3  4  5  4  3
         2  1  2  3  4  3  3  0  3  4  5  4
         3  2  1  2  3  4  4  3  0  3  4  5
         4  3  2  1  2  3  5  4  3  0  3  4
         3  4  3  2  1  2  4  5  4  3  0  3
         2  3  4  3  2  1  3  4  5  4  3  0]
   
    @test benzene.D == D
    
end