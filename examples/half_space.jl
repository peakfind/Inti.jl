using Inti
using LinearAlgebra
using StaticArrays
using SpecialFunctions: hankelh1

function helmholtz_kernel(target, source, k)
    x, y = Inti.coords(target), Inti.coords(source)
    yc = SVector(y[1], -y[2])
    d, dc = norm(x - y), norm(x - yc)
    d == 0 ? zero(ComplexF64) : im / 4*(hankelh1(0, k * d) - hankelh1(0, k * dc))
end

k = 50π
λ = 2π/k 
meshsize = λ/10

geo = Inti.parametric_curve(s -> SVector(cos(s), 2 + sin(s)), 0, 2π)
Γ = Inti.Domain(geo)
msh = Inti.meshgen(Γ; meshsize)
Q = Inti.Quadrature(msh; qorder = 5)

# Use a new kernel
K = let k = k
    (t, q) -> helmholtz_kernel(t, q, k) 
end
Sop = Inti.IntegralOperator(K, Q, Q)

# Compression
S₀ = Inti.assemble_matrix(Sop)

# Correction by using adaptive_correction
δS = Inti.adaptive_correction(Sop; tol = 1e-4, maxdist = 5*meshsize)

# Add δS to S₀
axpy!(1.0, δS, S₀)