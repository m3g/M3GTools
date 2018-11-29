#
# Function that wraps the coordinates to obtain minimum images
# around a defined center
# 
# It modifies the x,y,z input vectors
#

function wrap!(sides,x,y,z;center=[0.,0.,0.],sel=[-1])

  # If no selection was set, wrap all

  if sel[1] == -1 
    n = length(x)
    for i in 1:n
      wrapone!(sides,x[i],y[i],z[i],center)
    end


  # Wrap only the selection
  else
    n = length(sel)
    for i in 1:n
      wrapone!(sides,x[sel[i]],y[sel[i]],z[sel[i]],center)
    end
  end

end

function wrapone!(sides,x,y,z,center)

  x = x - center[1]
  y = y - center[2]
  z = z - center[3]

  x = x%sides[1]
  y = y%sides[2]
  z = z%sides[3]

  if x > sides[1]/2 ; x = x - sides[1] ; end
  if y > sides[2]/2 ; y = y - sides[2] ; end
  if z > sides[3]/2 ; z = z - sides[3] ; end

  if x < -sides[1]/2 ; x = x + sides[1] ; end
  if y < -sides[2]/2 ; y = y + sides[2] ; end
  if z < -sides[3]/2 ; z = z + sides[3] ; end

  x = x + center[1]
  y = y + center[2]
  z = z + center[3]

end
