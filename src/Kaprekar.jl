function kaprekar(a::Int;c::Int=0)
	a_str = if a<10
		"000$a"
	elseif 10<=a<100
		"00$a"
	elseif  100<=a<1000
		"0$a"
	else
		"$a"
	end
	a_arr = [parse(Int,s) for s in a_str]
	a_unq = unique(a_arr)
	min_len = length(a_unq)
	@assert min_len>=2 "Number must not be repdigits"
	n = length(a_arr)
	@assert n==4 "Number must four-digit integer counting leading zeros"
	kaprekar(a_arr,c)
end

function kaprekar(a::Vector{Int},c::Int)
	if a==[6,1,7,4] 
		return a, c
	end
	a_asc = sort(a)
	a_des = a_asc[end:-1:1]
	a_n = parse(Int,join(a_des)) - parse(Int,join(a_asc))
	kaprekar(a_n;c=c+1)
end

