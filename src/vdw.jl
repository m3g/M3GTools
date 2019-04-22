#
# This functions compute VdW interactions between two selections
#

#
# Computes the electrostatic interaction of a pair of atoms, given the distance
#

# Given the parameters for each atom

function vdwpair(d2 :: Float64, eps1 :: Float64, rmin1 :: Float64, 
                                eps2 :: Float64, rmin2 :: Float64)
  epspair = sqrt( eps1*eps2 )
  rminpair6 = (rmin1 + rmin2)^6
  return vdwpair( d2, epspair, rminpair6 )
end

# If the pair parameters were computed already

function vdwpair(d2 :: Float64, epspair :: Float64, rminpair6 :: Float64 )
  d = sqrt(d2)
  p6 = rminpair6 / d2^3
  p12 = p6*p6
  vdwpair = epspair*( p12 - 2*p6 )
  return vdwpair
end

#
# For a trajectory: x, y, z coordinates read from a DCD are given
#

# First, the general function with all parameters

function vdw(atoms :: Vector{Atom}, sel1 :: Vector{Int}, sel2 :: Vector{Int}, 
             sides :: Vector{Float64}, 
             x :: Vector{Float32}, y :: Vector{Float32}, z :: Vector{Float32}, 
             usecutoff :: Bool, cutoff :: Float64, pbc :: Bool)
  n1 = length(sel1)
  n2 = length(sel2)
  cutoff2 = cutoff^2
  vdw = 0.
  for i in 1:n1
    center = [ x[sel1[i]], y[sel1[i]], z[sel1[i]] ]
    if pbc 
      wrap!(sides,x,y,z,center=center,sel=sel2)
    end
    for j in 1:n2
      d2 = (x[sel1[i]] - x[sel2[j]])^2 + (y[sel1[i]] - y[sel2[j]])^2 + (z[sel1[i]] - z[sel2[j]])^2
      epspair = sqrt( atoms[i].eps*atoms[j].eps )
      rminpair6 = ( atoms[i].rmin + atoms[j].rmin )^6
      vdw = vdw + vdwpair(d2,epspair,rminpair6)
      if usecutoff
        if d2 < cutoff2
          vdw = vdw - vdwpair(cutoff2,epspair,rminpair6)
        end
      end
    end
  end
  return coulomb
end

# For a lazy user which does not want to setup the cutoff and the pbc, use defaults 

vdw(atoms :: Vector{Atom}, sel1 :: Vector{Int}, sel2 :: Vector{Int}, 
    sides :: Vector{Float64}, 
    x :: Vector{Float32}, y :: Vector{Float32}, z :: Vector{Float32} ) =
  vdw(atoms, sel1, sel2, sides, x, y, z, true, 12., true)

# If the user does not provide the sides, than we will not use PBCs

vdw(atoms :: Vector{Atom}, sel1 :: Vector{Int}, sel2 :: Vector{Int}, 
        x :: Vector{Float32}, y :: Vector{Float32}, z :: Vector{Float32} ) =
  vdw(atoms, sel1, sel2, [0., 0., 0.], x, y, z, true, 12., false)

#
# For computing interaction energies from the coordinates available in the atoms array
#

function vdw2( atoms :: Vector{Atom}, sel1 :: Vector{Int}, sel2 :: Vector{Int}, 
               usecutoff :: Bool, cutoff :: Float64 )
  n1 = length(sel1)
  n2 = length(sel2)
  vdw = 0.
  cutoff2 = cutoff^2
  for isel1 in 1:n1
    i = sel1[isel1]
    for isel2 in 1:n2
      j = sel2[isel2]
      d2 = (atoms[i].coor[1] - atoms[j].coor[1])^2 +
           (atoms[i].coor[2] - atoms[j].coor[2])^2 +
           (atoms[i].coor[3] - atoms[j].coor[3])^2 
      epspair = sqrt( atoms[i].eps*atoms[j].eps )
      rminpair6 = ( atoms[i].rmin + atoms[j].rmin )^6
      vdw = vdw + vdwpair(d2,epspair,rminpair6)
      if usecutoff
        if d2 < cutoff2
          vdw = vdw - vdwpair(cutoff2,epspair,rminpair6)
        end
      end
    end
  end
  return vdw
end

# Mapping to use a single function name

vdw(  atoms :: Vector{Atom}, sel1 :: Vector{Int}, sel2 :: Vector{Int}, usecutoff :: Bool, cutoff :: Float64 ) = 
  vdw2( atoms, sel1, sel2, usecutoff, cutoff )

# If the cutoff is not set, use default 

vdw( atoms :: Vector{Atom}, sel1 :: Vector{Int}, sel2 :: Vector{Int} ) = 
  vdw2( atoms, sel1, sel2, true, 12. )

