#
# Structures required for the computation of the 3D density of the gmd function
#

struct GMD3DGrid
  x :: Vector{Float64}
  atom :: Int64
  dmin :: Float64
  rho :: Float64
end

mutable struct MutableGMD3DGrid
  x :: Vector{Float64}
  atom :: Int64
  dmin :: Float64
  rho :: Float64
end

MutableGMD3DGrid() = MutableGMD3DGrid(zeros(3),0,0.,0.)

