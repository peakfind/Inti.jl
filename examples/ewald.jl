using SpecialFunctions: expint, erfc
using StaticArrays

"""
    ewald(x, y, k, α, d, a, N, M, J)

Compute the quasi-periodic Green functions by Ewald's method

# Arguments

- `k`: the wave number
- `α`: the quasi-momentum
- `d`: the period
- `a`: the parameter introduced due to contour deformation
- `N`, `M`, `J`: the number of truncated terms
"""
function ewald(x, y, k, α, d, a, N, M, J)
    spec = 0.0 + 0.0im
    spat = 0.0 + 0.0im
    # Precompute the reusable constants
    p = 2π/d
    q = a/d
    coef = (k*d/(2a))^2 # will be used in the spatial part

    # Spectrum part
    for n in -N:N
        # Precompute the reusable constants
        αₙ = α + n*p
        # γₙ = αₙ^2 >= k^2 ? sqrt(αₙ^2 - k^2) : -im*sqrt(k^2 - αₙ^2)
        γₙ = -im*sqrt(Complex(k^2 - αₙ^2))
        carg = γₙ*d/(2a)
        arg1 = carg + q*y
        arg2 = carg - q*y

        spec += exp(im*αₙ*x)*(exp(γₙ*y)*erfc(arg1) + exp(-γₙ*y)*erfc(arg2))/γₙ
    end

    # Spatial part
    for m in -M:M
        # Precompute the reusable constants
        xₘ = x - m*d
        rₘ = hypot(xₘ, y)
        arg = (rₘ*q)^2
        term = 0.0 + 0.0im

        for j in 0:J
            term += (coef^j)*expint(j + 1, arg)/factorial(j)
        end
        spat += term*exp(im*α*m*d)
    end

    return -spec/(4d) - spat/(4π)
end

"""
    gradx_ewald(x, y, k, α, d, a, N, M, J)

Compute the gradient of quasi-periodic Green functions by Ewald's method

# Arguments

- `k`: the wave number
- `α`: the quasi-momentum
- `d`: the period
- `a`: the parameter introduced due to contour deformation
- `N`, `M`, `J`: the number of truncated terms
"""
function gradx_ewald(x, y, k, α, d, a, N, M, J)
    ∂₁Gspec = 0.0 + 0.0im
    ∂₂Gspec = 0.0 + 0.0im
    ∂₁Gspat = 0.0 + 0.0im
    ∂₂Gspat = 0.0 + 0.0im
    # Precompute the reusable constants
    p = 2π/d
    q = a/d
    q′ = 2*q/sqrt(π) # used in ∂₂Gspec
    coef = (k*d/(2a))^2 # used in the spatial part

    # Spectrum part
    for n in -N:N
        # Precompute the reusable constants
        αₙ = α + n*p
        γₙ = -im*sqrt(Complex(k^2 - αₙ^2))
        carg = γₙ*d/(2a)
        arg₁ = carg + q*y
        arg₂ = carg - q*y
        arg₁² = arg₁^2
        arg₂² = arg₂^2

        # Reusable factors
        αₑ = exp(im*αₙ*x)
        γₑ⁺ = exp(γₙ*y)
        γₑ⁻ = exp(-γₙ*y)

        term1 = γₑ⁺*erfc(arg₁)
        term2 = γₑ⁻*erfc(arg₂)
        term3 = exp(γₙ*y)*exp(-arg₁²)
        term4 = exp(-γₙ*y)*exp(-arg₂²)

        ∂₁Gspec += im*αₙ*αₑ*(term1 + term2)/γₙ
        ∂₂Gspec += (αₑ*(term1 - term2) - q′*αₑ*(term3 - term4)/γₙ)
    end

    # Spatial part
    for m in -M:M
        # Precompute the reusable constants
        xₘ = x - m*d
        rₘ = hypot(xₘ, y)
        arg = (rₘ*q)^2
        term = 0.0 + 0.0im

        for j in 0:J
            term += (coef^j)*expint(j, arg)/factorial(j)
        end

        ∂₁Gspat += term*xₘ*exp(im*α*m*d)
        ∂₂Gspat += term*exp(im*α*m*d)
    end

    ∂₁G = (q^2)*∂₁Gspat/(2π) - ∂₁Gspec/(4d)
    ∂₂G = y*(q^2)*∂₂Gspat/(2π) - ∂₂Gspec/(4d)
    
    return SVector{2,ComplexF64}(∂₁G, ∂₂G)
end

#= quasi_periodic_helmoltz(x, y) = _quasi_periodic_helmholtz(x, y, :ewaldp)
periodic_helmhotz(x, y) = _periodic_helmholtz(x, y, :ewaldp) =#

#= function _quasi_periodic_helmholtz(x, y, method::Symbol)
    if method == :ewald
        return ewald()
    elseif method == :ewaldp
        return ewald_p()
    end
end =#

#= function _periodic_helmholtz(x, y, method::Symbol)
    if method == :ewald
        return ewald()
    elseif method == :ewaldp
        return ewald_p()
    end
end =#

#function ewald_p()
#end



