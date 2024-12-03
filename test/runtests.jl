using SynergisticGaussianSystems
using Test
using Distributions

# import here functions to test that are not exported in SynergisticGaussianSystems.jl
import SynergisticGaussianSystems: _bitstring_to_bitvec, _bitvec_to_bitstring, _vec_to_cov, _cov_to_vec, _cov_to_bitvec

# test important functions
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

    # test multivariate N = 4 process (synthetic data generated using MvNormal)
    local Σ = _vec_to_cov([.4,.4,.4,-.2,-.1,-.4]);
    local X = rand(MvNormal(zeros(size(Σ)[1]), Σ), 1000000);
    # true values computed analytically in Mathematica
    @test entropy_gaussian(X) - 4.258947420779911 < 10^(-1)
    @test total_correlation_gaussian(X) - 1.416806712038779 < 10^(-1)
    @test O_information_gaussian(X) - (-1.736742694823721) < 10^(-1)
end

# test helper functions
@testset "utils.jl" begin
    @test _bitstring_to_bitvec("001") == [0,0,1]
    @test _bitvec_to_bitstring([0,0,1]) == "001"
    local mat = [
        1.0  0.2  0.2;
        0.2  1.0  0.2;
        0.2  0.2  1.0;
        ]
    @test _vec_to_cov([.2,.2,.2]) == mat
    @test _cov_to_vec(mat) == [.2,.2,.2]
    @test _cov_to_bitvec([
		1 .3 .3;
		.3 1 -.2;
		.3 -.2 1
	]) == [0, 0, 1]
end

# use the following from the REPL to test the whole package
# (SynergisticGaussianSystems) pkg> test SynergisticGaussianSystems

# TODO:
# synch with TravisCI, appveyor, coverall etc so they stop saying it doesn't work