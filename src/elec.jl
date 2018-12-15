#
# This function computes electrostatic interaction between two selections
#

include("./qpair.jl")

function elec(sel1,sel2,charge,sides,x,y,z;cutoff=-1.,pbc=true)

  n1 = length(sel1)
  n2 = length(sel2)

  elec = 0.
  for i in 1:n1
    center = [ x[sel1[i]], y[sel1[i]], z[sel1[i]] ]
    if pbc 
      wrap!(sides,x,y,z,center=center,sel=sel2)
    end
    for j in 1:n2
      d = sqrt((x[sel1[i]] - x[sel2[j]])^2 + (y[sel1[i]] - y[sel2[j]])^2 + (z[sel1[i]] - z[sel2[j]])^2)
      if cutoff > 0.
        if d < cutoff
          elec = elec + qpair(d,charge[sel1[i]],charge[sel2[j]])
          shift = qpair(cutoff,charge[sel1[i]],charge[sel2[j]])
          elec = elec - shift
        end
      else
        elec = elec + qpair(d,charge[sel1[i]],charge[sel2[j]])
      end
    end
  end
  
  return elec
end

