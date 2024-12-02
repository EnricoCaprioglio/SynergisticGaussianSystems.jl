using SynergisticGaussianSystems
using Test

# import here functions to test that are not exported in SynergisticGaussianSystems.jl
import SynergisticGaussianSystems._bitstring_to_bitvec

@testset "SynergisticGaussianSystems.jl" begin
    # Write your tests here.
    @test _bitstring_to_bitvec("001") == [0,0,1]
    @test _bitstring_to_bitvec("011") == [0,1,1]
end

# use the following from the REPL to test the whole package
# (SynergisticGaussianSystems) pkg> test SynergisticGaussianSystems.jl