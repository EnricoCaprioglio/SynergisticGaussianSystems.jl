using LinearAlgebra
using Distributions

"""
Compute Gaussian Entropy from covariance matrix
"""
function _entropy_from_covariance(Σ::AbstractMatrix)
    @assert size(Σ)[1] == size(Σ)[2] "Σ is not a square matrix"
    @assert isposdef(Σ) "Σ is not a positive definite matrix"
    N = size(Σ)[1]
   return (1/2) * log(det(Σ) * (2*π*ℯ)^N) 
end

"""
Compute Gaussian Entropy from variance
"""
function _entropy_from_variance(σ::Real)::Real
    return (1/2) * log(σ * (2*π*ℯ))
end

"""
	function entropy_gaussian(X::AbstractArray; bias = false)

`X` can either be:

- your data, where each row corresponds to the data of a single process. Then `X` has to be a matrix with `N = # of rows` and `no_datapoints = # of columns`

- a covariance matrix, then `size(X)[1] == size(X)[2]`

- data of a single proces, then `X` is a vector

- if `X` is a number, then that number has to be the variance of a Gaussian process
"""
function entropy_gaussian(X; bias = false)::Real
	
	if length(size(X)) > 1 # then X has multiple processes
		if size(X)[1] == size(X)[2] # then this is a covariance matrix
			
            Σ = X
			return _entropy_from_covariance(Σ)
			
		else # then X is a mutlivariate Gaussian system
			
			# compute covariance matrix along the columns
			Σ = cov(X, dims = 2)
			return _entropy_from_covariance(Σ)
		end
		
	elseif length(X) > 1 # then X is a sigle process timeseries
		
        σ = var(X)
		return _entropy_from_variance(σ)

    elseif 0 ≤ X ≤ 1 # then X is a variance
		
        return _entropy_from_variance(X)

	end
	
end


"""
	function function total_correlation(X::AbstractArray)

`X` can either be:

- your data, where each row corresponds to the data of a single process. Then `X` has to be a matrix with `N = # of rows` and `no_datapoints = # of columns`

- a covariance matrix, then `X` has to satisfy `size(X)[1] == size(X)[2]`
"""
function total_correlation_gaussian(X::AbstractArray)::Real

	N = size(X)[1]

	if size(X)[1] == size(X)[2] # then X is a covariance Σ
		Σ = X
		# get individual processes entropy
		individual_H = [entropy_gaussian(Σ[i,i]) for i in 1:N]
		# get entropy of the whole system
		H_whole = entropy_gaussian(Σ)

		return sum(individual_H) - H_whole
	else
		
		# get individual processes entropy
		individual_H = [entropy_gaussian(X[i, :]) for i in 1:N]
		# get entropy of the whole system
		H_whole = entropy_gaussian(X)
	
		return sum(individual_H) - H_whole
	end
end

"""
    function O_information_gaussian(X::AbstractArray)

`X` can either be:

- your data, where each row corresponds to the data of a single process. Then `X` has to be a matrix with `N = # of rows` and `no_datapoints = # of columns`

- a covariance matrix, then `X` has to satisfy `size(X)[1] == size(X)[2]`
"""
function O_information_gaussian(X::AbstractArray)::Real

	N = size(X)[1]
	residuals = zeros(N)

	if size(X)[1] == size(X)[2]

		Σ = X
		
		for i in 1:N
			# get covariance of subsystem X₋ᵢ and compute entropy
			Σ_subsys = Σ[setdiff(1:end, i), setdiff(1:end, i)]
			# compute total correlation of subsystem X₋ᵢ
			residuals[i] = total_correlation_gaussian(Σ_subsys)
		end

		# compute TC of the whole system
		TC = total_correlation_gaussian(Σ)

		# compute the O-information Ω = (∑ᵢ TC(X₋ᵢ)) - (N-2) TC(X)
		return sum(residuals) - N*TC + 2*TC

	else

		for i in 1:N
			# isolate subsystem X₋ᵢ
			X_minus_i = X[setdiff(1:end, i), :]
			# compute total correlation of subsystem X₋ᵢ
			residuals[i] = total_correlation_gaussian(X_minus_i)
		end

		# compute TC of the whole system
		TC = total_correlation_gaussian(X)

		# compute the O-information Ω = (∑ᵢ TC(X₋ᵢ)) - (N-2) TC(X)
		return sum(residuals) - (N -2) * TC
		
	end
end