#=
# qpgreen.jl
# 
# Plot the real part of the quasi-periodic Green functions
# by Ewald's method
# 
# Reference: Mathematical and Computational Methods in Photonics 
#            and Phononics, Fig 2.11
=#

using CairoMakie
include("ewald.jl")

function field_real(x, y)
    k = 2 + 0.5im # wave number
    α = π/8       # quasi-momentum
    a = sqrt(π)   # parameter introduced by Ewald's method
    d = 1         # period
    N = 5
    M = 5
    J = 10
    return real(ewald(x, y, k, α, d, a, N, M, J))
end

xs = range(-5, 5, length=200)
ys = range(-5, 5, length=200)
zs = [field_real(x, y) for x in xs, y in ys]

fig, ax, hm = heatmap(xs, ys, zs, colormap=:winter)
Colorbar(fig[:, end+1], hm)
fig