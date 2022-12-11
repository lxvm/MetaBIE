"""
    potential

Compute the quasi-periodic Green's function for the Helmholtz equation using a
window function
"""
function potential(Δx, Δy, k, kx, L, m, x0=1e-3*m*L, x1=m*L)
    Δϕ = cis(-kx*L)
    ϕn = Δϕ^0
    acc = window(Δx,x0,x1)*ϕn*hankelh1(0, k*sqrt((Δx+0*L)^2 + Δy^2))
    for n in 1:m
        ϕn *= Δϕ
        acc += window(Δx+n*L,x0,x1)*hankelh1(0, k*sqrt((Δx+n*L)^2 + Δy^2))*ϕn
            +  window(Δx-n*L,x0,x1)*hankelh1(0, k*sqrt((Δx-n*L)^2 + Δy^2))*ϕn'
    end
    return 0.25im*acc
end

function window(x::T, x0::T, x1::T) where {T<:AbstractFloat}
    absx = abs(x)
    if absx <= x0
        return T(1)
    elseif absx >= x1
        return T(0)
    else
        u = (absx - x0)/(x1 - x0)
        return exp(T(2)*exp(-inv(u))/(u-T(1)))
    end
end

function single_layer_op(k, kx, L, Γ, targets, m)
    S = Matrix{ComplexF64}(undef, length(targets), length(Γ))
    for (j,r2) in enumerate(Γ), (i,r1) in enumerate(targets)
        Δr = r1-r2
        ϕ = potential(Δr[1], Δr[2], k, kx, L, m)
        isfinite(ϕ) || @show r1 r2
        S[i,j] = ϕ
    end
    return S
end