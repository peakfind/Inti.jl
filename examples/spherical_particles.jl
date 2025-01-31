#=
# spherical_particles.jl
# 
# Use Inti.jl to solve a scattering problem of a periodic 
# array of spherical particles 
# 
# Reference: Mathematical and Computational Methods in Photonics 
#            and Phononics, Section 4.6
=#

using GLMakie

# TODO: implement struct Medium
ω = 1
ϵ₁ = 1
μ₁ = 1
ϵ₂ = 5
μ₂ = 1
θ = pi/8

k₁ = ω*sqrt(ϵ₁*μ₁)
k₂ = ω*sqrt(ϵ₂*μ₂)

α = k₁*sin(θ)
β = k₁*cos(θ)

xs = LinRange(-3*pi, 3*pi, 100)
ys = LinRange(-3*pi, 3*pi, 100)
u = [3*cos(α*x - β*y) for x in xs, y in ys]

surface(x, y, u)
