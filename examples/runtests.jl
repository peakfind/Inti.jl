using Test
include("ewald.jl")

"""
Tests for the function ewald(x, y, k, α, d, a, N, M, J)

# Parameters for Tables 2 - 5

## Table 2: x = 0.01*d, y = 0, k = 2/d, α = sqrt(2)/d, d = 2π
- test 1: a = 2, N = 3, M = 2, J = 7 −0.4595298795 - 0.3509130869im
- test 2: a = 2, N = 2, M = 1, J = 4 −0.4595297462 - 0.3509130866im

## Table 3: x = 0.01*d, y = 0, k = 10/d, α = 5*sqrt(2)/d, d = 2π
- test 1:  a = 1  N = 3  M = 2 J = 78 −0.3538185358 − 0.1769332429im
- test 2:  a = 2  N = 4  M = 2 J = 28 −0.3538172307 − 0.1769332383im
- test 3:  a = 3  N = 5  M = 1 J = 18 −0.3538172307 − 0.1769332383im
- test 4:  a = 4  N = 6  M = 1 J = 13 −0.3538172307 − 0.1769332383im
- test 5:  a = 5  N = 7  M = 0 J = 11 −0.3538172307 − 0.1769332383im
- test 6:  a = 6  N = 9  M = 0 J = 9  −0.3538172307 − 0.1769332383im
- test 7:  a = 7  N = 10 M = 0 J = 8  −0.3538172307 − 0.1769332383im
- test 8:  a = 8  N = 11 M = 0 J = 8  −0.3538172307 − 0.1769332383im
- test 9:  a = 9  N = 12 M = 0 J = 7  −0.3538172307 − 0.1769332383im
- test 10: a = 10 N = 13 M = 0 J = 7  −0.3538172307 − 0.1769332383im

## Table 4: x = 0.01*d, y = 0, k = 2/d, α = 3/d, d = 2π
- test 1: a = 2 N = 3 M = 2 J = 7 −0.7617463954 - 0.0006387129177im

## Table 5: x = 0.5*d, k = 2/d, α = sqrt(2)/d, d = 2π
- test 1: y = 0.1*d a = 2 N = 2 M = 2 J = 7 0.3306805081 − 0.1778394385im
- test 2: y = 0.5*d a = 2 N = 3 M = 2 J = 7 0.3596087433 − 0.04626396800im
"""

@testset verbose = true "Numerical results for Ewald's method in Linton1998" begin
    @testset "Tests for Table 2" begin
        d = 2π
        x = 0.01*d
        y = 0
        k = 2/d
        α = sqrt(2)/d
        @test ewald(x, y, k, α, d, 2, 3, 2, 7) ≈ −0.4595298795 - 0.3509130869im atol=1e-10
        @test ewald(x, y, k, α, d, 2, 2, 1, 4) ≈ −0.4595297462 - 0.3509130866im atol=1e-10
        
    end
    @testset "Tests for Table 3" begin
        d = 2π
        x = 0.01*d
        y = 0
        k = 10/d
        α = 5*sqrt(2)/d
        @test ewald(x, y, k, α, d, 1, 3, 2, 78) ≈ −0.3538185358 − 0.1769332429im atol=1e-10
        @test ewald(x, y, k, α, d, 2, 4, 2, 28) ≈ −0.3538172307 − 0.1769332383im atol=1e-10
        @test ewald(x, y, k, α, d, 3, 5, 1, 18) ≈ −0.3538172307 − 0.1769332383im atol=1e-10
        @test ewald(x, y, k, α, d, 4, 6, 1, 13) ≈ −0.3538172307 − 0.1769332383im atol=1e-10
        @test ewald(x, y, k, α, d, 5, 7, 0, 11) ≈ −0.3538172307 − 0.1769332383im atol=1e-10
        @test ewald(x, y, k, α, d, 6, 9, 0, 9) ≈ −0.3538172307 − 0.1769332383im atol=1e-10
        @test ewald(x, y, k, α, d, 7, 10, 0, 8) ≈ −0.3538172307 − 0.1769332383im atol=1e-10
        @test ewald(x, y, k, α, d, 8, 11, 0, 8) ≈ −0.3538172307 − 0.1769332383im atol=1e-10
        @test ewald(x, y, k, α, d, 9, 12, 0, 7) ≈ −0.3538172307 − 0.1769332383im atol=1e-10
        @test ewald(x, y, k, α, d, 10, 13, 0, 7) ≈ −0.3538172307 − 0.1769332383im atol=1e-10
    end
    @testset "Tests for Table 4" begin
        d = 2π
        x = 0.01*d
        y = 0
        k = 2/d
        α = 3/d
        @test ewald(x, y, k, α, d, 2, 3, 2, 7) ≈ −0.7617463954 - 0.0006387129177im atol=1e-10
    end
    @testset "Tests for Table 5" begin
        d = 2π 
        x = 0.5*d
        k = 2/d
        α = sqrt(2)/d
        @test ewald(x, 0.1*d, k, α, d, 2, 2, 2, 7) ≈ 0.3306805081 − 0.1778394385im atol=1e-10
        @test ewald(x, 0.5*d, k, α, d, 2, 3, 2, 7) ≈ 0.3596087433 − 0.04626396800im atol=1e-10
    end
end
nothing