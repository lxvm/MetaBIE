module MetaBIE

using LinearAlgebra

using StaticArrays
using Bessels
using LinearMaps
using IterativeSolvers
using FFTW

include("fourier.jl")
include("kernel.jl")
include("solve.jl")

end # module MetaBIE
