#
# Structures required for the computation of the 3D density
#

struct Density
  x :: Vector{Float64}
  atom :: Int64
  dmin :: Float64
  rho :: Float64
end

mutable struct MutableDensity
  x :: Vector{Float64}
  atom :: Int64
  dmin :: Float64
  rho :: Float64
end

MutableDensity() = MutableDensity(zeros(3),0,0.,0.)

