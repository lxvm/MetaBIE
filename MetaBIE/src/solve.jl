struct MetaProblem{TA,TB}
    L::Float64 # period
    k::Float64 # incident wavenumber
    θ::Float64 # incident angle
    α::TA # metasurface parameter
    β::TB # metasurface parameter
    n::Int # Nystrom points
    m::Int # window points
end

struct MetaSolution{T}
    L::Float64 # period
    k::Float64 # incident wavenumber
    θ::Float64 # incident angle
    σ1::Vector{T} # charge density in half-space 1
    σ2::Vector{T} # charge density in half-space 2
    n::Int # Nystrom points
    m::Int # window points
end

function directsolve(prob::MetaProblem)
    Γx = range(0.0, 2pi, length=2prob.n+1)[1:2prob.n]
    # Γ = map(x -> SVector{2,Float64}(x, zero(x)), Γx)
    k = prob.k
    kz, kx = k .* sincospi(prob.θ)
    S = single_layer_op(k, kx, prob.L, prob.n, prob.m)
    αs = map(x -> prob.α.(prob.L*x/2pi), Γx)
    βs = map(x -> prob.β.(prob.L*x/2pi), Γx)
    ui = map(x -> cispi(kx*prob.L*x/2pi), Γx)
    ξ1 = (-0.5I + Diagonal(im*k*αs)*S)\(-2im*k*Diagonal(αs)*ui)
    ξ2 = (-0.5I + Diagonal(im*k*βs)*S)\(-2im*kz*ui)
    MetaSolution(prob.L, prob.k, prob.θ, 0.5(ξ1+ξ2), 0.5(ξ1-ξ2), prob.n, prob.m)
end

function fastsolve(prob::MetaProblem)
    Γx = range(0.0, 2pi, length=2prob.n+1)[1:2prob.n]
    # Γ = map(x -> SVector{2,Float64}(x, zero(x)), Γx)
    k = prob.k
    kz, kx = k .* sincospi(prob.θ)
    DFT = plan_fft!(rand(ComplexF64, 2prob.n))
    ŝ = DFT*single_layer_op_circ(k, kx, prob.L, prob.n, prob.m)
    αs = map(x -> prob.α.(prob.L*x/2pi), Γx)
    βs = map(x -> prob.β.(prob.L*x/2pi), Γx)
    ui = map(x -> cispi(kx*prob.L*x/2pi), Γx)
    f1 = let k=k, αs=αs, ŝ=ŝ
        x -> -0.5x + DFT*(Diagonal(im*k*αs)*(DFT\(Diagonal(ŝ)*x)))
    end
    ξ1 = DFT\gmres(LinearMap(f1, 2prob.n), DFT*(-2im*k*Diagonal(αs)*ui))
    f2 = let k=k, βs=βs, ŝ=ŝ
        x -> -0.5x + DFT*(Diagonal(im*k*βs)*(DFT\(Diagonal(ŝ)*x)))
    end
    ξ2 = DFT\gmres(LinearMap(f2, 2prob.n), DFT*(-2im*kz*ui))
    MetaSolution(prob.L, prob.k, prob.θ, 0.5(ξ1+ξ2), 0.5(ξ1-ξ2), prob.n, prob.m)
end



function (sol::MetaSolution)(x, z)
    z == 0 && return ComplexF64(NaN, NaN)
    k = sol.k
    kz, kx = k .* sincospi(sol.θ)
    Γx = range(-sol.L/2, sol.L/2, length=2sol.n+1)[1:2sol.n]
    Γ = map(x -> SVector{2,Float64}(x, zero(x)), Γx)
    S = single_layer_op(k, kx, sol.L, Γ, Ref(SVector(x,z)), sol.m)
    ui = cis(kx*x+kz*z)
    if z > 0
        # @show ui S * sol.σ1
        return ui + only(S * sol.σ1)
    else
        # @show ui S * sol.σ2
        return ui + only(S * sol.σ2)
    end
end

function iswoodanomaly(L,kx,k)
    Δk = 2pi/L
    kxnp = kx + zero(Δk)
    kxnm = kx - zero(Δk)
    while true
        kxnp^2 > k^2 && kxnm^2 > k^2 && return false
        kxnp^2 ≈ k^2 && return true
        kxnp += Δk
        kxnm^2 ≈ k^2 && return true
        kxnm -= Δk
    end
end