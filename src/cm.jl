#
# This function computes the center of mass of a selection
#

function cm(selection,mass,x,y,z)

  cm = [ 0., 0., 0. ]
  totmass = 0.
  for i in 1:length(selection)
    cm[1] = cm[1] + mass[selection[i]]*x[selection[i]] 
    cm[2] = cm[2] + mass[selection[i]]*y[selection[i]] 
    cm[3] = cm[3] + mass[selection[i]]*z[selection[i]] 
    totmass = totmass + mass[selection[i]]
  end
  cm[1] = cm[1] / totmass
  cm[2] = cm[2] / totmass
  cm[3] = cm[3] / totmass

  return cm
end

