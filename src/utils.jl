function _bitstring_to_bitvec(bin_string)
	bitvec = []
	for s in bin_string
		push!(bitvec, parse(Int, s))
	end
	return bitvec
end