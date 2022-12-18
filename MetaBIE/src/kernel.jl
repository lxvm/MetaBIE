"""
    potential

Compute the quasi-periodic Green's function for the Helmholtz equation using a
window function
"""
function potential(Δx, Δy, k, kx, L, m, x0=0.5*m*L, x1=m*L)
    Δϕ = cis(-kx*L)
    ϕn = Δϕ^0
    acc = window(Δx,x0,x1)*ϕn*hankelh1(0, k*sqrt((Δx+0*L)^2 + Δy^2))
    for n in 1:m
        ϕn *= Δϕ
        acc +=  window(Δx+n*L,x0,x1)*hankelh1(0, k*sqrt((Δx+n*L)^2 + Δy^2))*ϕn +
                window(Δx-n*L,x0,x1)*hankelh1(0, k*sqrt((Δx-n*L)^2 + Δy^2))*ϕn'
    end
    return 0.25im*acc
end

"Window function introduced in https://doi.org/10.1016/j.jcp.2013.12.047"
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

const nystrom_cache = Dict{Int,Vector{Float64}}()

function cachedrule(n)
    haskey(nystrom_cache, n) && return nystrom_cache[n]
    nystrom_cache[n] = nystrom_weights(n)
end

"Nyström weights on PTR grid"
function nystrom_weights(n)
    R = Vector{Float64}(undef, 2n)
    for i in 1:2n
        r = cospi(i-1)/2n
        for m in 1:n-1
            r += cospi(m*(i-1)/n)/m
        end
        R[i] = -2pi*r/n
    end
    return R
end

"Nystrom method & Martensen-Kussmaul quadrature on ptr grid"
function single_layer_op!(S, k, kx, L, N, m, jmax=2N, x0=1e-3*2pi*m, x1=2pi*m)
    Δs = L/2pi # Jacobian of rescaling [0,L] -> [0,2π]
    k_ = Δs*k
    kx_ = Δs*kx
    Δx = pi/N # step size in PTR grid

    R = cachedrule(N)
    Δϕ = cispi(-2kx_)

    for j in 1:jmax, i in 1:2N
        Δij = if -N <= (dij = i-j) <= N
            dij
        elseif dij > 0
            dij - 2N
        else
            2N + dij
        end
        Δxij = Δx*Δij
        kr = k_*abs(Δxij)
        
        ϕn = Δϕ^0
        g = zero(Δϕ)
        for n in 1:m
            ϕn *= Δϕ
            g += window(Δxij+n*2pi,x0,x1)*hankelh1(0, k_*abs(Δxij+n*2pi))*ϕn +
                 window(Δxij-n*2pi,x0,x1)*hankelh1(0, k_*abs(Δxij-n*2pi))*ϕn'
        end
        g *= 0.25im

        # apply MK kernel for the singular n=0 term in G_QP sum
        s1 = -besselj0(kr)/4pi
        s2 = if i == j
            (0.25im - γ/2pi - log(0.5k_)/2pi)
        else
            0.25im*window(Δxij,x0,x1)*hankelh1(0,kr) - s1*log(4*sin(0.5Δxij)^2)
        end
        S[i,j] = Δs*(R[abs(Δij)+1]*s1 + Δx*(s2 + g))
    end
end

const γ = 0.57721566490153286060651209008240243104215933593992 # Euler-Mascheroni constant

function single_layer_op(k, kx, L, n, m)
    S = Matrix{ComplexF64}(undef, 2n, 2n)
    if isinteger(L*kx/2pi)
        single_layer_op!(S, k, kx, L, n, m, 1)
        for i in 2:2n
            circshift!(view(S, :, i), view(S, :, 1), i-1)
        end
    else
        single_layer_op!(S, k, kx, L, n, m)
    end
    return S
end

function single_layer_op_circ(k, kx, L, n, m)
    @assert isinteger(L*kx/2pi) "Cannot use fast method in quasiperiodic case"
    S = Vector{ComplexF64}(undef, 2n)
    single_layer_op!(S, k, kx, L, n, m, 1)
    return S
end

function single_layer_op(k, kx, L, Γ, targets, m)
    S = Matrix{ComplexF64}(undef, length(targets), length(Γ))
    Δs = L/length(Γ)
    for (j,r2) in enumerate(Γ), (i, r1) in enumerate(targets)
        dx, Δy = r1-r2
        Δx = if -0.5L <= dx <= 0.5L
            dx
        elseif dx > 0
            dx - L
        else
            L + dx
        end
        S[i,j] = Δs*potential(Δx, Δy, k, kx, L, m)
    end
    return S
end