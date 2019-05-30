#
# This function computes the center of mass of a selection
#

function cm(selection,simulation,x,y,z)

  cm = [ 0., 0., 0. ]
  totmass = 0.
  for i in 1:length(selection)
    cm[1] = cm[1] + simulation.atoms[selection[i]].mass*x[selection[i]] 
    cm[2] = cm[2] + simulation.atoms[selection[i]].mass*y[selection[i]] 
    cm[3] = cm[3] + simulation.atoms[selection[i]].mass*z[selection[i]] 
    totmass = totmass + simulation.atoms[selection[i]].mass
  end
  cm[1] = cm[1] / totmass
  cm[2] = cm[2] / totmass
  cm[3] = cm[3] / totmass

  return cm
end

