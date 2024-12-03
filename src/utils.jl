using LinearAlgebra
using Distributions
using Random

function _bitstring_to_bitvec(bin_string)
	bitvec = []
	for s in bin_string
		push!(bitvec, parse(Int, s))
	end
	return bitvec
end

function _bitvec_to_bitstring(bitvec)
	bit_string = ""
	for v in bitvec
		bit_string = bit_string*string(v)
	end
	return bit_string
end

function _vec_to_cov(v::AbstractArray)

    # find desired matrix size n = N*(N-1)/2
    Ns = collect(2:11); ns = map(x-> Int(x*(x-1)/2), Ns)
	n = length(v); n_id = findfirst(x-> x == n, ns)
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