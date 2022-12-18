include("params.jl")
include("test.jl")

# explicitly set the number of threads to 1
BLAS.set_num_threads(1)
FFTW.set_num_threads(1)

plt = plot_prob(prob)
savefig(plt, "params.png")

sol = MetaBIE.directsolve(prob)

plt = plot_sol(sol)
savefig(plt, "densities.png")

plt = plot_sol_field(sol, 2)
savefig(plt, "fields.png")

Ns = 2 .^ (1:12) # numbers of grid points to use
M  = 20 # number of terms to include in windowed sum of quasiperiodic Green's function

dsols = check_convergence(Ns, MetaBIE.directsolve, M)
fsols = check_convergence(Ns, MetaBIE.fastsolve, M)

plt = convergence_plot(dsols, Ns)
plot!(plt, log2.(Ns), Ns .^ (-4); color=3, ls=:dash, label="O(N⁻⁴)", ylim=(1e-15,1))
savefig(plt, "densities_conv_direct.png")

plt = field_convergence_plot(dsols, Ns)
plot!(plt, log2.(Ns), Ns .^ (-4); color=3, ls=:dash, label="O(N⁻⁴)", ylim=(1e-15,1))
savefig(plt, "fields_conv_direct.png")

plt = convergence_plot(fsols, Ns)
plot!(plt, log2.(Ns), Ns .^ (-4); color=3, ls=:dash, label="O(N⁻⁴)", ylim=(1e-15,1))
savefig(plt, "densities_conv_fast.png")

plt = field_convergence_plot(fsols, Ns)
plot!(plt, log2.(Ns), Ns .^ (-4); color=3, ls=:dash, label="O(N⁻⁴)", ylim=(1e-15,1))
savefig(plt, "fields_conv_fast.png")

dtimes = check_scaling(Ns, MetaBIE.directsolve, M)
ftimes = check_scaling(Ns, MetaBIE.fastsolve, M)

plt = plot(; xguide="N", yguide="wall time (s)", scale=:log10, legend=:topleft, title="Method timing comparison")
plot!(plt, Ns, dtimes; label="direct", color=1, markershape=:x)
plot!(Ns, N -> 1e-8N^3; color=1, ls=:dash, label="O(N³)")
plot!(plt, Ns, ftimes; label="fast", color=2, markershape=:x)
plot!(Ns, N -> 5e-6N*log(N); color=2, ls=:dash, label="O(Nlog(N))")
savefig(plt, "scaling.png")
