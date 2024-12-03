using SynergisticGaussianSystems
using Test
using Distributions

# import here functions to test that are not exported in SynergisticGaussianSystems.jl
import SynergisticGaussianSystems._bitstring_to_bitvec
import SynergisticGaussianSystems._bitvec_to_bitstring
import SynergisticGaussianSystems._vec_to_cov

@testset "Gaussian_IT.jl" begin
    
    # test single process
    @test isapprox(entropy_gaussian(randn(1000000)), (1/2) * log((2*π*ℯ)); atol=0.01)

    # test N = 3 example
    local Σ = _vec_to_cov([.4,.4,-.2]);
    # true values computed analytically in Mathematica
    @test entropy_gaussian(Σ) - 3.980991790470895 < 10^(-12)
    @test total_correlation_gaussian(Σ) - 0.2758238091431229 < 10^(-12)
    @test O_information_gaussian(Σ) - (-0.08105942473821748) < 10^(-12)

    # test N = 4 example
    local Σ = _vec_to_cov([.4,.4,.4,-.2,-.1,-.4]);
    # true values computed analytically in Mathematica
    @test entropy_gaussian(Σ) - 4.258947420779911 < 10^(-12)
    @test total_correlation_gaussian(Σ) - 1.416806712038779 < 10^(-12)
    @test O_information_gaussian(Σ) - (-1.736742694823721) < 10^(-12)

    # test multivariate N = 4 process
    local Σ = _vec_to_cov([.4,.4,.4,-.2,-.1,-.4]);
    local X = rand(MvNormal(zeros(size(Σ)[1]), Σ), 100000);
    # true values computed analytically in Mathematica
    @test entropy_gaussian(X) - 4.258947420779911 < 10^(-2)
    @test total_correlation_gaussian(X) - 1.416806712038779 < 10^(-2)
    @test O_information_gaussian(X) - (-1.736742694823721) < 10^(-2)
end

# use the following from the REPL to test the whole package
# (SynergisticGaussianSystems) pkg> test SynergisticGaussianSystems

# TODO:
# synch with TravisCI, appveyor, coverall etc so they stop saying it doesn't work