using MetaBIE


L = 2pi # domain period
k = 12.3 # wavenumber
θi = -pi/2 # incidence angle (measured wrt x axis)
θt = pi/3 # desired angle of transmission

n = 100 # number of terms in Nystrom grid
m = 20 # number of terms in windowed sum

kx = k*cos(θi) # parallel wavenumber
kz = k*sin(θi) # perpendicular wavenumber
@assert kz != 0 "Wood anomaly"
d = cos(θt) - kz # metasurface parameter
c = -sin(θi) # metasurface parameter
# L = 2pi/(k*d) # period of designer metasurface

α = MetaBIE.FourierSeries([0, c, c], L)
β = MetaBIE.FourierSeries([0, c,-c], L)

prob = MetaBIE.MetaProblem(L, k, θi, α, β, n, m)