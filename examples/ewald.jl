using SpecialFunctions: expint, erfc

"""
    ewald(x, y, k, α, d, a, N, M, J)

Compute the quasi-periodic Green functions by Ewald's method

TODO: 
     1. Numerical Stability: sqrt ?
     2. Precomputation: Cache factorials, powers, and phase factorials
     3. Maybe parallelize the code

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
        γₙ = αₙ^2 >= k^2 ? sqrt(αₙ^2 - k^2) : -im*sqrt(k^2 - αₙ^2)
        carg = γₙ*d/(2a)
        arg1 = carg + q*y
        arg2 = carg - q*y

        spec += exp(im*αₙ*x)*(exp(γₙ*y)*erfc(arg1) + exp(-γₙ*y)*erfc(arg2))/γₙ
    end

    # spatial part
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



