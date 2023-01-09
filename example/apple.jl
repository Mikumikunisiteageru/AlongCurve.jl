# example/apple.jl

using AlongCurve
using PyPlot

image = imread(joinpath(pkgdir(AlongCurve), "example", "input.png"));

matrix = image[:, :, 1] .>= 0.5;

m, n = size(matrix)

@time coords = order(matrix) # 1.2 s

palette = [(1.0, 0.0, 0.0), (1.0, 0.5, 0.0), (1.0, 1.0, 0.0), (0.5, 1.0, 0.0), 
		   (0.0, 1.0, 0.0), (0.0, 1.0, 0.5), (0.0, 1.0, 1.0), (0.0, 0.5, 1.0), 
		   (0.0, 0.0, 1.0), (0.5, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 0.0, 0.5)]

colormat = ones(m, n, 3)
for i = eachindex(coords)
	colormat[coords[i].I..., 1:3] .= palette[mod1(i, 12)]
end

imsave(joinpath(pkgdir(AlongCurve), "example", "output.png"), colormat)
