using MetaBIE
using FourierSeriesEvaluators


# L = 2pi # domain period
k = 10.0 # wavenumber
θi = -1/2 # incidence angle (measured wrt x axis) (in fractions of radians)
θt = -1/6 # desired angle of transmission (in fractions of radians)

n = max(20, round(Int, 2pi/k)) # number of terms in Nystrom grid
m = 30 # number of terms in windowed sum

kx = k*cospi(θi) # parallel wavenumber
kz = k*sinpi(θi) # perpendicular wavenumber
# @assert kz != 0 "Wood anomaly"
d = cospi(θt) - cospi(θi) # metasurface parameter
c = -sinpi(θi) # metasurface parameter
L = 2pi/(k*abs(d)) # period of designer metasurface
@assert !MetaBIE.iswoodanomaly(L,kx,k) "Wood anomaly"

α = FourierSeries([c, c]; period=L, offset=-1)
β = FourierSeries([c,-c]; period=L, offset=-1)

# α = FourierSeries([0]; period=L, offset=-1)
# β = FourierSeries([0]; period=L, offset=-1)

prob = MetaBIE.MetaProblem(L, k, θi, α, β, n, m)