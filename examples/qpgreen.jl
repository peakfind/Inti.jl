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

function field_imag(x, y)
    k = 2 + 0.5im # wave number
    α = π/8       # quasi-momentum
    a = sqrt(π)   # parameter introduced by Ewald's method
    d = 1         # period
    N = 5
    M = 5
    J = 10
    return imag(ewald(x, y, k, α, d, a, N, M, J))
end

function grad1_field_real(x, y)
    k = 2 + 0.5im # wave number
    α = π/8       # quasi-momentum
    a = sqrt(π)   # parameter introduced by Ewald's method
    d = 1         # period
    N = 5
    M = 5
    J = 10
    ∇G = gradx_ewald(x, y, k, α, d, a, N, M, J)

    return real(∇G[1])
end

function grad2_field_real(x, y)
    k = 2 + 0.5im # wave number
    α = π/8       # quasi-momentum
    a = sqrt(π)   # parameter introduced by Ewald's method
    d = 1         # period
    N = 5
    M = 5
    J = 10
    ∇G = gradx_ewald(x, y, k, α, d, a, N, M, J)

    return real(∇G[2])
end

function grad1_field_imag(x, y)
    k = 2 + 0.5im # wave number
    α = π/8       # quasi-momentum
    a = sqrt(π)   # parameter introduced by Ewald's method
    d = 1         # period
    N = 5
    M = 5
    J = 10
    ∇G = gradx_ewald(x, y, k, α, d, a, N, M, J)

    return imag(∇G[1])
end

function grad2_field_imag(x, y)
    k = 2 + 0.5im # wave number
    α = π/8       # quasi-momentum
    a = sqrt(π)   # parameter introduced by Ewald's method
    d = 1         # period
    N = 5
    M = 5
    J = 10
    ∇G = gradx_ewald(x, y, k, α, d, a, N, M, J)

    return imag(∇G[2])
end

xs = range(-5, 5, length=200)
ys = range(-5, 5, length=200)
ℜG = [field_real(x, y) for x in xs, y in ys]
ℑG = [field_imag(x, y) for x in xs, y in ys]
ℜ∂G1 = [grad1_field_real(x, y) for x in xs, y in ys]
ℜ∂G2 = [grad2_field_real(x, y) for x in xs, y in ys]
ℑ∂G1 = [grad1_field_imag(x, y) for x in xs, y in ys]
ℑ∂G2 = [grad2_field_imag(x, y) for x in xs, y in ys]

fig = Figure()
ax₁₁ = Axis(fig[1, 1], title = L"$\Re G$")
ax₁₂ = Axis(fig[1, 2], title = L"$\Im G$")
ax₂₁ = Axis(fig[2, 1], title = L"\Re\left(\partial_{1}G\right)")
ax₂₂ = Axis(fig[2, 2], title = L"\Re\left(\partial_{2}G\right)")
ax₃₁ = Axis(fig[3, 1], title = L"\Im\left(\partial_{1}G\right)")
ax₃₂ = Axis(fig[3, 2], title = L"\Im\left(\partial_{2}G\right)")
heatmap!(ax₁₁, xs, ys, ℜG, colormap=:winter)
heatmap!(ax₁₂, xs, ys, ℑG, colormap=:winter)
heatmap!(ax₂₁, xs, ys, ℜ∂G1, colormap=:winter)
heatmap!(ax₂₂, xs, ys, ℜ∂G2, colormap=:winter)
heatmap!(ax₃₁, xs, ys, ℑ∂G1, colormap=:winter)
heatmap!(ax₃₂, xs, ys, ℑ∂G2, colormap=:winter)
fig