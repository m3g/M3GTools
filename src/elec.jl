#
# This function computes electrostatic interaction between two selections
#

function elec(sel1,sel2,charge,sides,x,y,z)

  Coulomb = 332.05382e0

  n1 = length(sel1)
  n2 = length(sel2)

  elec = 0.
  for i in 1:n1
    center = [ x[sel1[i]], y[sel1[i]], z[sel1[i]] ]
    wrap!(sides,x,y,z,center=center,sel=sel2)
    for j in 1:n
      d = sqrt((x[sel1[i]] - x[sel2[j]])^2 + (y[sel1[i]] - y[sel2[j]])^2 + (z[sel1[i]] - z[sel2[j]])^2)
      elec = elec + qpair(d,q1,q2)
    end
  end
  
  return elec
end

# Computes the electrostatic interaction of a pair of atoms, given the distance

function qpair(d,q1,q2)

  Coulomb = 332.05382e0
  qpair = Coulomb*q1*q2 / d

  return qpair
end
