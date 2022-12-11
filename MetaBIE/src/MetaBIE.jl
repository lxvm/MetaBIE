module MetaBIE

using LinearAlgebra

using StaticArrays
using Bessels

include("fourier.jl")
include("kernel.jl")
include("solve.jl")

end # module MetaBIE
