# test/runtests.jl

using AlongCurve
using Test

matrix = BitMatrix([0 1 1; 0 1 0; 1 0 0])
result = showorder(matrix)
@test result == [0 3 4; 0 2 0; 1 0 0] || result == [0 2 1; 0 3 0; 4 0 0]

@testset "Case 1" begin
	matrix = BitMatrix([
		0 1 1 0
		1 0 0 0 
		0 1 1 1
		0 0 1 0
		0 1 0 0])
	@test pathlength(matrix) == 3 + 4 * sqrt(2)
end

@testset "Case 2" begin
	matrix = BitMatrix([
		0 1 1 0 0
		1 0 0 0 1
		0 1 1 1 0
		0 0 1 0 0])
	@test pathlength(matrix) == 3 + 4 * sqrt(2)
end

@testset "Case 3" begin
	matrix = BitMatrix([
		0 0 1 0
		1 1 1 1
		0 0 1 0
		0 0 1 0])
	@test pathlength(matrix) == 4 + 2 * sqrt(2)
end

@testset "Case 4" begin
	matrix = BitMatrix([
		0 1 1
		0 1 0
		1 0 0])
	result = showorder(matrix)
	@test result == [0 3 4; 0 2 0; 1 0 0] || result == [0 2 1; 0 3 0; 4 0 0]
	@test pathlength(matrix) == 2 + 1 * sqrt(2)
end

@testset "Case 5" begin
	matrix = BitMatrix([
		0 1 1 0 0
		1 0 0 0 1
		0 1 1 1 0
		0 0 1 0 0
		0 1 0 0 0])
	@test_throws ErrorException order(matrix)
end
