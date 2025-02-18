using Test
include("muller.jl")

@testset verbose = true "Example in Mathematical and Computational Methods in Photonics and Phononics" begin
    function f(z)
        z = complex(z)
        return sin(z) + 5 + im
    end
    @test muller(f, -1.2, -1.4, -1.5) ≈ -1.369601247093987 - 2.313220941769529im atol=1e-10
end

@testset verbose = true "Examples in Roots.jl" begin
    @test muller(x->x^3-1, 0, 0.5im, 0.5) ≈ -0.5 + 0.8660254037im atol=1e-10
    @test muller(x->x^3-1, 0.5, 0.5im, -0.5) ≈ 1.0 atol=1e-10
    @test muller(x->x^2+2, 0.0, 0.5, 1.0) ≈ 0 - 1.4142135623im atol=1e-10
    @test muller(x->(x-5)*x*(x+5), 0.39, 0.65, 0.6) ≈ 0.00 atol=1e-10
end