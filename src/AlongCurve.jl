# src/AlongCurve.jl

module AlongCurve

export order, showorder, pathlength

using JuMP
using GLPK

const IJ = CartesianIndex{2}

offsets(matrix::BitMatrix, bij::IJ, sij::Vector{<:IJ}) = 
	filter(ij -> matrix[ij], 
		filter(ij -> checkbounds(Bool, matrix, ij), [bij] .+ sij))

adjacents(matrix::BitMatrix, ij::IJ) = offsets(matrix, ij, 
	[IJ(Int.(sign.(sincosd(k)))) for k = 0:45:315])

function simplify(matrix::BitMatrix)
	pixels = findall(matrix)
	n = length(pixels)
	dict = Dict(pixels .=> 1:n)
	adjres = [getindex.((dict,), a) for a = adjacents.((matrix,), pixels)]
	ends = findall(isone.(length.(adjres)))
	segments = Vector{Int}[]
	selected = falses(n)
	while ! isempty(ends)
		to = [pop!(ends)]
		segment = Int[]
		while isone(length(to))
			i = first(to)
			push!(segment, i)
			selected[i] = true
			to = filter(k -> ! selected[k], adjres[i])
		end
		push!(segments, segment)
	end
	return pixels, n, adjres, segments, selected
end

function dist(ij1::IJ, ij2::IJ, M::Float64)
	h = hypot((ij1-ij2).I...)
	return h > 1.5 ? M : h
end

function subsets(n::Int, d::Matrix{Float64})
	a = Vector{Int}[Int[]]
	for i = 1:n
		for s = a
			all(isfinite.(d[s, i])) && push!(a, vcat(s, i))
		end
	end
	return filter(s -> 2 <= length(s) <= n-1, a)
end

function hamilton(d::Matrix{Float64}; optimizer=GLPK)
	n = size(d, 1)
	model = direct_model(optimizer.Optimizer())
	@variable(model, x[1:n, 1:n], binary=true)
	@constraint(model, sum(x, dims=1) .== 1)
	@constraint(model, sum(x, dims=2) .== 1)
	for s = subsets(n, d)
		@constraint(model, sum(x[s, s]) <= length(s) - 1)
	end
	dbound = min.(d, n + 3.0)
	@objective(model, Min, sum(x .* dbound))
	optimize!(model)
	next = only.(findall.(eachcol(Bool.(value.(x)))))
	i = 1
	while (d[i, next[i]]) <= 1.5
		i = next[i]
	end
	path = Int[]
	for k = 1:n
		i = next[i]
		k < n && d[i, next[i]] > 1.5 && error("No solution!")
		push!(path, i)
	end
	return path
end

function order(matrix::BitMatrix; optimizer=GLPK)
	pixels, n, adjres, segments, selected = simplify(matrix)
	v0 = findall(.!selected)
	nv = length(v0)
	ns = length(segments)
	nt = nv + ns * 2
	v = copyto!(Vector{Int}(undef, nt), v0)
	d = Matrix{Float64}(undef, nt, nt)
	insegments = zeros(Int, n)
	for k = 1 : ns
		segment = segments[k]
		v[nv + 2*k .- [1, 0]] .= a, b = first(segment), last(segment)
		insegments[a] = insegments[b] = k
	end
	for i = 1 : nt
		d[i, i] = Inf
		for j = i+1 : nt
			d[i, j] = d[j, i] = dist(pixels[v[i]], pixels[v[j]], Inf)
		end
	end
	for i = nv+1 : 2 : nt
		d[i, i+1] = d[i+1, i] = 0.0
	end
	ham = hamilton(d; optimizer=optimizer)
	path = Int[]
	i = 0
	while i < nt
		i += 1
		h = v[ham[i]]
		k = insegments[h]
		if k > 0
			if h == first(segments[k])
				@assert v[ham[i+1]] == last(segments[k])
				append!(path, segments[k])
			else # h == last(segments[k])
				@assert v[ham[i+1]] == first(segments[k])
				append!(path, reverse(segments[k]))
			end
			i += 1
		else
			push!(path, h)
		end
	end
	return pixels[path]
end

function pathlength(matrix::BitMatrix; optimizer=GLPK)
	path = order(matrix; optimizer=optimizer)
	return sum(dist(path[i], path[i+1], Inf) for i = 1:length(path)-1)
end
	
function showorder(matrix::BitMatrix; optimizer=GLPK)
	path = order(matrix; optimizer=optimizer)
	m, n = size(matrix)
	newmat = zeros(Int, m, n)
	newmat[path] .= eachindex(path)
	return newmat
end

end # module AlongCurve
