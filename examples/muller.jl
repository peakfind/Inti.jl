"""
    muller(f, xₙ₋₂,  xₙ₋₁, xₙ; xtol, ftol, maxit)

Finding a zero of a function defined on the complex plane

# Arguments

- `f`: nonlinear function
- `xₙ₋₂`, `xₙ₋₁`, `xₙ`: three initial approximations
- `xtol`: tolerance of sucessive approxiamte roots
- `ftol`: tolerance of absolute value of the function f at approxiamte root
- `maxit`: maximum number of iterations

# References

- Numerical Methods in Scientific Computing Vol. I, Section 6.2.3
- [Roots.jl](https://github.com/JuliaMath/Roots.jl)
"""
function muller(f::Function, xₙ₋₂,  xₙ₋₁, xₙ; xtol=1e-5, ftol=1e-12, maxit=100)
    # Promote initial points to complex for type stability
    xₙ₋₂, xₙ₋₁, xₙ = complex(xₙ₋₂), complex(xₙ₋₁), complex(xₙ)

    # The initial approximations should be distinct
    if xₙ₋₂ == xₙ₋₁ || xₙ₋₁ == xₙ || xₙ₋₂ == xₙ
        throw(ArgumentError("xₙ₋₂, xₙ₋₁, xₙ should be distinct!"))
    end

    fₙ₋₁ = f(xₙ₋₁)
    fₙ₋₂ = f(xₙ₋₂)

    for iter in 1:maxit
        fₙ = f(xₙ)

        # Precompute reusable constants
        q = (xₙ - xₙ₋₁)/(xₙ₋₁ - xₙ₋₂)
        q² = q^2
        q1 = q + one(q)

        A = q*fₙ - q*q1*fₙ₋₁ + q²*fₙ₋₂
        B = (q + q1)*fₙ - q1^2*fₙ₋₁ + q²*fₙ₋₂
        C = q1 * fₙ
        
        # Compute the discriminant
        Δ = B^2 - 4*A*C

        # Determine the denominator
        deno⁺ = B + sqrt(Δ)
        deno⁻ = B - sqrt(Δ)
        deno = abs(deno⁺) > abs(deno⁻) ? deno⁺ : deno⁻

        x = xₙ - 2*(xₙ - xₙ₋₁)*fₙ/deno

        # Update the approximations
        xₙ₋₂, xₙ₋₁, xₙ  = xₙ₋₁, xₙ, x
        fₙ₋₂, fₙ₋₁= fₙ₋₁, fₙ

        # Check the termination criteria
        fx = f(x)
        if abs(xₙ - xₙ₋₁) < xtol && abs(fx) < ftol
            println("Number of iterations: $iter")
            return xₙ
        end
    end

    @warn "Maximum iterations have been used. The solution may not converge."
    return xₙ 
end