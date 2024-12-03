using LinearAlgebra
using Distributions
using Random

"""
    function _bitstring_to_bitvec(bin_string)

Given a string of zeros and ones it outputs a vector of zeros and ones.

e.g.,

	_bitstring_to_bitvec("001")
	→ [0, 0, 1]
"""
function _bitstring_to_bitvec(bin_string)
	bitvec = []
	for s in bin_string
		push!(bitvec, parse(Int, s))
	end
	return bitvec
end

"""
    function _bitvec_to_bitstring(bitvec)
        
Given a vector of zeros and ones it outputs a string of zeros and ones.

e.g.,

	_bitvec_to_bitstring([0,0,1])
	→ "001"
"""
function _bitvec_to_bitstring(bitvec)
	bit_string = ""
	for v in bitvec
		bit_string = bit_string*string(v)
	end
	return bit_string
end

"""
	_vec_to_cov(vec::AbstractArray)

Construct covariance matrix from vector `vec` with elements between `-1` and `1`.
This function does not check whether the covariance is positive definite or not.

e.g.,

	_vec_to_cov([.2,.2,.2])
	→ 3×3 Matrix{Float64}:
		 1.0  0.2  0.2
		 0.2  1.0  0.2
		 0.2  0.2  1.0
"""
function _vec_to_cov(v::AbstractArray)

    @assert all(x -> abs(x) >= 0 && abs(x) <= 1, v) "All elements must have absolute values between 0 and 1"
    # find desired matrix size n = N*(N-1)/2
    Ns = collect(2:11); ns = map(x -> Int(x*(x-1)/2), Ns)
	n = length(v); n_id = findfirst(x -> x == n, ns)
    @assert !isnothing(n_id) "Could not construct Σ, the length of v is likely incorrect or too large"
	N = Ns[n_id]
	
    # construct Σ
	Σ = zeros(N, N); id = 1
	for i in 1:N
		Σ[i,i] = 1
		for j in i+1:N
			Σ[i,j] = v[id]
			Σ[j,i] = Σ[i,j]
			id += 1
		end
	end

	return Σ
end

"""
	_cov_to_vec(mat::AbstractArray)

Output vector with off-diagonal elements of `mat::AbstractArray`.

e.g.,

	_cov_to_vec([
		1 1 1;
		1 1 0;
		1 0 1
	])
	→ [1, 1, 0]
"""
function _cov_to_vec(mat::AbstractArray)
	return [mat[i,j] for i in 1:size(mat, 1) for j in i+1:size(mat, 2)]
end

"""
	_cov_to_bitvec(mat::AbstractArray)
Output configuration string of covariance matrix `mat::AbstractArray`.
This function uses the configuration convection `i > 0 ? 0 : 1`. Do not confuse with `_cov_to_vec`.

e.g.,

	_cov_to_bitvec([
		1 .3 .3;
		.3 1 -.2;
		.3 -.2 1
	])
	→ [0, 0, 1]
"""
function _cov_to_bitvec(mat::AbstractArray)
	vec = _cov_to_vec(mat)
	bitvec = [i > 0 ? 0 : 1 for i in vec]
	return bitvec
end