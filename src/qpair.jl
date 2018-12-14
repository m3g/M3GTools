#
# Computes the electrostatic interaction of a pair of atoms, given the distance
#

function qpair(d,q1,q2)

  Coulomb = 332.05382e0
  qpair = Coulomb*q1*q2 / d

  return qpair
end
