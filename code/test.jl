using LinearAlgebra
using Random

using StaticArrays
using Plots
using FFTW
using BenchmarkTools
using FourierSeriesEvaluators

using MetaBIE

function check_convergence(Ns, solver, m=20)
    k = 10.0 # wavenumber
    θi = -1/2 # incidence angle (measured wrt x axis)
    θt = -1/6 # desired angle of transmission

    kx = k*cospi(θi) # parallel wavenumber
    kz = k*sinpi(θi) # perpendicular wavenumber
    d = cospi(θt) - cospi(θi) # metasurface parameter
    c = -sinpi(θi) # metasurface parameter
    L = 2pi/(k*abs(d)) # period of designer metasurface
    @assert !MetaBIE.iswoodanomaly(L,kx,k) "Wood anomaly"

    α = FourierSeries([c, c]; period=L, offset=-1)
    β = FourierSeries([c,-c]; period=L, offset=-1)

    sols = Vector{MetaBIE.MetaSolution{ComplexF64}}(undef, length(Ns))
    @info "starting convergence run for $(nameof(solver))"
    for (i,N) in enumerate(Ns)
        prob = MetaBIE.MetaProblem(L, k, θi, α, β, N, m)
        @info "starting $N"
        t = time()
        sols[i] = solver(prob)
        @info "finished $N in $(time()-t) (s)"
    end
    @info "finished convergence run"
    return sols
end

function convergence_plot(sols_, Ns_)
    sols = view(sols_, 1:length(sols_)-1)
    Ns = view(Ns_, 1:length(Ns_)-1)
    refsol = sols_[end]
    refbase = Int(log2(Ns_[end]))
    bases = [Int(log2(N)) for N in Ns]
    errs1 = [maximum(abs, sols[i].σ1 - refsol.σ1[1:(2^(refbase-bases[i])):end]) for i in eachindex(sols)]
    errs2 = [maximum(abs, sols[i].σ2 - refsol.σ2[1:(2^(refbase-bases[i])):end]) for i in eachindex(sols)]
    plot(; xguide="log2(N)", yguide="L∞ error relative to ||σ||∞", yscale=:log10, title="Self-convergence of surface densities")
    plot!(bases, errs1 ./ maximum(abs, refsol.σ1); label="σ1", markershape=:x)
    plot!(bases, errs2 ./ maximum(abs, refsol.σ2); label="σ2", markershape=:x)
end

function field_convergence_plot(sols, Ns, y=1.0, p=501)
    y_ = y*2pi/sols[1].k
    xs = range(0.0, sols[1].L, length=p)
    tsol1 = map(x -> sols[end](x, y_), xs)
    tsol2 = map(x -> sols[end](x, -y_), xs)
    errs1 = map(sol -> maximum(abs, tsol1 .- map(x -> sol(x, y_), xs)), first(sols, length(sols)-1))
    errs2 = map(sol -> maximum(abs, tsol2 .- map(x -> sol(x, -y_), xs)), first(sols, length(sols)-1))
    plot(; xguide="log2(N)", yguide="L∞ error relative to ||u||∞", yscale=:log10, title="Self-convergence of total fields")
    plot!(Int.(log2.(first(Ns, length(Ns)-1))), errs1 ./ maximum(abs, tsol1); label="u1", markershape=:x)
    plot!(Int.(log2.(first(Ns, length(Ns)-1))), errs2 ./ maximum(abs, tsol2); label="u2", markershape=:x)
end

function check_scaling(Ns, solver, m=20)
    k = 10.0 # wavenumber
    θi = -1/2 # incidence angle (measured wrt x axis)
    θt = -1/6 # desired angle of transmission

    kx = k*cospi(θi) # parallel wavenumber
    kz = k*sinpi(θi) # perpendicular wavenumber
    d = cospi(θt) - cospi(θi) # metasurface parameter
    c = -sinpi(θi) # metasurface parameter
    L = 2pi/(k*abs(d)) # period of designer metasurface
    @assert !MetaBIE.iswoodanomaly(L,kx,k) "Wood anomaly"

    α = FourierSeries([c, c]; period=L, offset=-1)
    β = FourierSeries([c,-c]; period=L, offset=-1)

    times = Vector{Float64}(undef, length(Ns))
    @info "starting timing run for $(nameof(solver))"
    for (i,N) in enumerate(Ns)
        prob = MetaBIE.MetaProblem(L, k, θi, α, β, N, m)
        @info "starting $N"
        t = time()
        times[i] = @belapsed $solver($prob)
        @info "finished $N in $(time()-t) (s)"
    end
    @info "finished timing run"
    return times
end

function plot_prob(prob)
    xs = range(0, prob.L, length=20prob.n)
    plot(; xguide="x", yguide="Metasurface parameter", xticks=([0, prob.L], ["0" "2π"]))
    plot!(xs, x -> real(prob.α(x)); label="re(α)", ls=:solid)
    plot!(xs, x -> imag(prob.α(x)); label="im(α)", ls=:dash)
    plot!(xs, x -> real(prob.β(x)); label="re(β)", ls=:solid)
    plot!(xs, x -> imag(prob.β(x)); label="im(β)", ls=:dash)
end

function plot_sol(sol)
    xs = range(0, sol.L, length=2sol.n+1)
    plot(; xguide="x", yguide="Surface density", xticks=([0, sol.L], ["0" "2π"]))
    σ1 = vcat(sol.σ1, sol.σ1[1])
    σ2 = vcat(sol.σ2, sol.σ2[1])
    plot!(xs, real.(σ1); label="re(σ1)", ls=:solid)
    plot!(xs, imag.(σ1); label="im(σ1)", ls=:dash)
    plot!(xs, real.(σ2); label="re(σ2)", ls=:solid)
    plot!(xs, imag.(σ2); label="im(σ2)", ls=:dash)
end

function plot_sol_field(sol, np=1, nλ=2, nx=200, nz=201)
    xs = range(-np*sol.L, np*sol.L, length=nx)
    zs = range(-nλ*2pi/sol.k, nλ*2pi/sol.k, length=nz)
    heatmap(xs, zs, (x,y) -> real(sol(x,y));
    xguide="x", yguide="z",
    xticks=([-np*sol.L, 0, np*sol.L], ["-$(2*np)L" "0" "$(2*np)L"]),
    yticks=([-nλ*2pi/sol.k, 0, nλ*2pi/sol.k], [ "-$(nλ)λ", "0", "$(nλ)λ"]),
    title="Re(u) field plot", aspectratio=1, color=:RdBu)
end

function show_rank(ks=[10.0, 100.0, 1000.0], kx=0.0, L=2pi, n=20, m=20, tol=1e-6, M=20, seed=18336)
    Random.seed!(seed)
    sources = map(x -> SVector(x, zero(x)), rand(n)*L/m)
    targets = map(x -> SVector(x, zero(x)), rand(n)*L/m)
    ranks = Vector{Int}(undef, m)
    plt = plot(;xticks=([0, L], ["0", "2π"]), xguide="x", yguide="relative rank")
    for k in ks
        for ii in 1:m
            ξ = SVector(L*(ii-1)/(m-1), 0.0)
            S = MetaBIE.single_layer_op(k, kx, L, sources .+ [ξ], targets, M)
            ranks[ii] = findfirst(<=(tol), svdvals(S))
        end
        plot!(plt, [L*(ii-1)/(m-1) for ii in 1:m], ranks; label="k=$k")
    end
    return plt
end