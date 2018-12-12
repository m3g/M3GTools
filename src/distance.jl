#
# This function computes distance between two selections
# (that is, the minimum distance between their atoms)
#

function distance(sel1,sel2,sides,x,y,z)

  n1 = length(sel1)
  n2 = length(sel2)

  distance = Inf
  for i in 1:n1
    center = [ x[sel1[i]], y[sel1[i]], z[sel1[i]] ]
    wrap!(sides,x,y,z,center=center,sel=sel2)
    for j in 1:n2
      distance = min(distance, (x[sel1[i]] - x[sel2[j]])^2 +
                               (y[sel1[i]] - y[sel2[j]])^2 +
                               (z[sel1[i]] - z[sel2[j]])^2 )
    end
  end 
  distance = sqrt(distance)

  return distance
end

