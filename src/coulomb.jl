#
# This functions compute electrostatic interaction between two selections
#

#
# Computes the electrostatic interaction of a pair of atoms, given the distance
#

function qpair(d,q1,q2)
  Coulomb = 332.05382e0
  qpair = Coulomb*q1*q2 / d
  return qpair
end

#
# For a trajectory: x, y, z coordinates read from a DCD are given
#

# First, the general function with all parameters

function coulomb(atoms :: Vector{Atom}, sel1 :: Vector{Int}, sel2 :: Vector{Int}, 
                 sides :: Vector{Float64}, 
                 x :: Vector{Float32}, y :: Vector{Float32}, z :: Vector{Float32}, 
                 usecutoff :: Bool, cutoff :: Float64, pbc :: Bool)
  n1 = length(sel1)
  n2 = length(sel2)
  coulomb = 0.
  for i in 1:n1
    center = [ x[sel1[i]], y[sel1[i]], z[sel1[i]] ]
    if pbc 
      wrap!(sides,x,y,z,center=center,sel=sel2)
    end
    for j in 1:n2
      d = sqrt((x[sel1[i]] - x[sel2[j]])^2 + (y[sel1[i]] - y[sel2[j]])^2 + (z[sel1[i]] - z[sel2[j]])^2)
      coulomb = coulomb + qpair(d,atoms[sel1[i]].charge,atoms[sel2[j]].charge)
      if usecutoff
        if d < cutoff
          coulomb = coulomb - qpair(cutoff,atoms[sel1[i]].charge,atoms[sel2[j]].charge)
        end
      end
    end
  end
  return coulomb
end

# For a lazy user which does not want to setup the cutoff and the pbc, use defaults 

coulomb(atoms :: Vector{Atom}, sel1 :: Vector{Int}, sel2 :: Vector{Int}, 
        sides :: Vector{Float64}, 
        x :: Vector{Float32}, y :: Vector{Float32}, z :: Vector{Float32} ) =
  coulomb(atoms, sel1, sel2, sides, x, y, z, true, 12., true)

# If the user does not provide the sides, than we will not use PBCs

coulomb(atoms :: Vector{Atom}, sel1 :: Vector{Int}, sel2 :: Vector{Int}, 
        x :: Vector{Float32}, y :: Vector{Float32}, z :: Vector{Float32} ) =
  coulomb(atoms, sel1, sel2, [0., 0., 0.], x, y, z, true, 12., false)

#
# For computing interaction energies from the coordinates available in the atoms array
#

function coulomb2( atoms :: Vector{Atom}, sel1 :: Vector{Int}, sel2 :: Vector{Int}, usecutoff :: Bool, cutoff :: Float64 )
  n1 = length(sel1)
  n2 = length(sel2)
  coulomb = 0.
  for isel1 in 1:n1
    i = sel1[isel1]
    for isel2 in 1:n2
      j = sel2[isel2]
      d = sqrt( (atoms[i].coor[1] - atoms[j].coor[1])^2 +
                (atoms[i].coor[2] - atoms[j].coor[2])^2 +
                (atoms[i].coor[3] - atoms[j].coor[3])^2 )
      coulomb = coulomb + qpair(d,atoms[i].charge,atoms[j].charge)
      if usecutoff 
        if d < cutoff
          coulomb = coulomb - qpair(cutoff,atoms[i].charge,atoms[j].charge)
        end
      end 
    end
  end
  return coulomb
end

# Mapping to use a single function name

coulomb(  atoms :: Vector{Atom}, sel1 :: Vector{Int}, sel2 :: Vector{Int}, usecutoff :: Bool, cutoff :: Float64 ) = 
  coulomb2( atoms, sel1, sel2, usecutoff, cutoff )

# If the cutoff is not set, use default 

coulomb( atoms :: Vector{Atom}, sel1 :: Vector{Int}, sel2 :: Vector{Int} ) = 
  coulomb( atoms, sel1, sel2, true, 12. )

