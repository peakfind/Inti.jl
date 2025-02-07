#=
# spherical_particles.jl
# 
# Use Inti.jl to solve a scattering problem of a periodic 
# array of spherical particles 
# 
# Reference: Mathematical and Computational Methods in Photonics 
#            and Phononics, Section 4.6
=#

#= TODO: some notions that will be used in Biebs.jl
1. fundamental solutions (kernels)
   laplacian, periodic_laplacian, quasi_periodic_laplacian
   helmholtz, periodic_helmholtz, quasi_periodic_helmholtz
2. layer operator
   single, double and the other two (confirm their terms)
   ∫G(x, y)φ(y) ds(y)
3. numerical methods of the integral equation (numie)
   nystrom method: nystrom.jl
   projection method: projection.jl

   numie/nystrom.jl, numie/projection.jl
4. we should use wrapers
   for example: 
   When we compute the kernel functions, we should consider two types:
   direct probelms and the eigenvalue problems
   
   we can use the same inner functions with prefix _ 
   _somefunc()
   and call different functions when we do different computations (julia's multiple dispatch)
   somefunc(x, y) (or target, source) # for direct problem
   somefunc(x, y, param)              # for eigenvalue problem

=#

#= 
# Visualize the scattered field
# xs = range(-3*π, 3*π, 100)
# ys = range(-3*π, 3*π, 100)
# uˢ = map(, Iterators.product(xs, ys))
# uⁱ
# uᵗ = uⁱ + uˢ
# fig, ax, hm = heatmap(xs, ys, , colormap = :inferno)
=#

using Inti
using StaticArrays
using LinearAlgebra
using GLMakie
include("ewald.jl")

struct Medium
   ϵ::Float64
   μ::Float64
end    

function get_wavenumber(ω, m::Medium)
   return ω*sqrt(m.ϵ*m.μ)
end   

Ω₁ = Medium(1, 1)
Ω₂ = Medium(5, 1)
ω = 1   # frequency
θ = π/8 # incident angle
k₁ = get_wavenumber(ω, Ω₁)
k₂ = get_wavenumber(ω, Ω₂)
α = k₁*sin(θ)
β = k₁*cos(θ)
d = 1 # period

# Create a mesh for the geometry
∂Ω₂ = Inti.parametric_curve(θ->SVector(0.4*cos(θ), 0.4*sin(θ)), 0.0, 2π, labels = ["∂Ω₂"])
msh = Inti.meshgen(∂Ω₂; meshsize = π/8)

# Create a quadrature
Q = Inti.Quadrature(msh; qorder = 3)

# Construct rhs 
rhs₁ = map(Q) do q
   x = q.coords
   return 3*exp(im*(α*x[1] - β*x[2]))
end

rhs₂ = map(Q) do q
    x = q.coords
    ν = q.normal
    return 3*exp(im*(α*x[1] - β*x[2]))*im*(α*ν[1] - β*ν[2])/Ω₁.μ
end

append!(rhs₁, rhs₂)

# Integral operators
# quasi-periodic Green functions
function quasi_periodic_helmholtz(target, source, k, α, d)
   t = Inti.coords(target)
   s = Inti.coords(source)
   dist = norm(t - s)
   # Parameters used by Ewald's method
   a = 
   N = 
   M = 
   J = 

   # the singularity at t = s needs to be handled separately, 
   # so just put a zero now
   if dist == 0
      return zero(ComplexF64) 
   else 
      return ewald() # TODO:
   end
end

# Single layer operators
G₁ = let k = k₁, α = α, d = d
   (t, s) -> quasi_periodic_helmholtz(t, s, k, α, d)
end
G₂ = let k = k₂
   (t, s) -> quasi_periodic_helmholtz(t, s, k)
end
S₁ = Inti.IntegralOperator(G₁, Q, Q)
S₂ = Inti.IntegralOperator(G₂, Q, Q)

# Adjoint double layer operators


A = []
# Solve the linear system
φ = A \ rhs₁