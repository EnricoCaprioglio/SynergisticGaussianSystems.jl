module SynergisticGaussianSystems

# Write your package code here.
include("utils.jl")

include("Gaussian_IT.jl")
export entropy_gaussian
export total_correlation_gaussian
export O_information_gaussian

end
