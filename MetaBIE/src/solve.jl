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

function solve(prob::MetaProblem)
    Γx = range(-prob.L/2, prob.L/2, length=prob.n)
    Γ = map(x -> SVector{2,Float64}(x, zero(x)), Γx)
    k = prob.k
    kz, kx = k .* sincos(prob.θ)
    S = single_layer_op(k, kx, prob.L, Γ .+ [SVector(step(Γx)/2, 0)], Γ, prob.m)
    αs = prob.α.(Γx)
    βs = prob.β.(Γx)
    ui = map(x -> cis(k*x), Γx)
    ξ1 = (-0.5I + Diagonal(im*k*αs)*S)\(2im*k*Diagonal(αs)*ui)
    ξ2 = (-0.5I + Diagonal(im*k*βs)*S)\(-2im*kz*ui)
    MetaSolution(prob.L, prob.k, prob.θ, 0.5(ξ1+ξ2), 0.5(ξ1-ξ2), prob.n, prob.m)
end