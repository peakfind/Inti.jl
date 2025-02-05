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

using GLMakie
using Inti
using StaticArrays
using Meshes
using LinearAlgebra

# Set medium parameters
# ω = 1
# ϵ₁ = 1
# μ₁ = 1 
# ϵ₂ = 5
# μ₂ = 1
# θ = pi/8

# Get the wavenumbers
# k₁ = ω*sqrt(ϵ₁*μ₁)
# k₂ = ω*sqrt(ϵ₂*μ₂)

# α = k₁*sin(θ)
# β = k₁*cos(θ)

# xs = LinRange(-3*pi, 3*pi, 100)
# ys = LinRange(-3*pi, 3*pi, 100)
# u = [3*cos(α*x - β*y) for x in xs, y in ys]

# surface(xs, ys, u)

# Create a mesh for the geometry
∂Ω = Inti.parametric_curve(θ->SVector(0.4*cos(θ), 0.4*sin(θ)), 0.0, 2*pi, labels = ["∂Ω"])
msh = Inti.meshgen(∂Ω; meshsize = 0.05)
viz(msh)

# Create a quadrature
Q = Inti.Quadrature(msh; qorder = 3)

# Assemble the Matrix and rhs vector
# using single_layer or something else to construct discrete integral operator

# Create the block operator

# solve the density function

# Visualize the scattered field
xs = range(-3*π, 3*π, 100)
ys = range(-3*π, 3*π, 100)
uˢ = map(, Iterators.product(xs, ys))
# uⁱ
# uᵗ = uⁱ + uˢ
fig, ax, hm = heatmap(xs, ys, , colormap = :inferno)

function quasi_periodic_helmoltz()
end