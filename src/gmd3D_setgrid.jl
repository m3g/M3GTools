#
# This function builds the grid of the 3D density function and fills an array of
# mutable structures of type MutableGMD3DGrid, containing the position of the atoms of 
# grid, the closest atom to that position, and distance. 
#

function gmd3D_setgrid(mysim,solute;dmax=7.0,step=0.5,dmin=1.0)

  xmin = [ mysim.atoms[solute[1]].coor[1], 
           mysim.atoms[solute[1]].coor[2], 
           mysim.atoms[solute[1]].coor[3] ] 
  xmax = [ mysim.atoms[solute[1]].coor[1], 
           mysim.atoms[solute[1]].coor[2], 
           mysim.atoms[solute[1]].coor[3] ] 

  for i in 2:length(solute)
    xmin[1] = min(xmin[1],mysim.atoms[solute[i]].coor[1])
    xmin[2] = min(xmin[2],mysim.atoms[solute[i]].coor[2])
    xmin[3] = min(xmin[3],mysim.atoms[solute[i]].coor[3])
    xmax[1] = max(xmax[1],mysim.atoms[solute[i]].coor[1])
    xmax[2] = max(xmax[2],mysim.atoms[solute[i]].coor[2])
    xmax[3] = max(xmax[3],mysim.atoms[solute[i]].coor[3])
  end

  xmin[1] = xmin[1] - dmax
  xmin[2] = xmin[2] - dmax
  xmin[3] = xmin[3] - dmax
  xmax[1] = xmax[1] + dmax
  xmax[2] = xmax[2] + dmax
  xmax[3] = xmax[3] + dmax

  nx = (xmax[1]-xmin[1])/step+1
  ny = (xmax[2]-xmin[2])/step+1
  nz = (xmax[3]-xmin[3])/step+1

  # Counting the number of grid points 

  dmin2 = dmin^2
  dmax2 = dmax^2
  ngrid = 0
  for ix in 1:nx
    x = xmin[1] + step*(ix-1)
    for iy in 1:ny
      y = xmin[2] + step*(iy-1)
      for iz in 1:nz
        z = xmin[3] + step*(iz-1)
        for iat in 1:length(solute)
          d2 = ( mysim.atoms[solute[iat]].coor[1] - x )^2 + 
               ( mysim.atoms[solute[iat]].coor[2] - y )^2 + 
               ( mysim.atoms[solute[iat]].coor[3] - z )^2  
          if d2 <= dmax2 && d2 > dmin2
            ngrid = ngrid + 1
            break
          end
        end
      end
    end
  end

  # Initializing grid structure

  grid = Vector{MutableGMD3DGrid}(undef,ngrid)
 
  # Saving the data for each grid point
      
  igrid = 0
  for ix in 1:nx
    x = xmin[1] + step*(ix-1)
    for iy in 1:ny
      y = xmin[2] + step*(iy-1)
      for iz in 1:nz
        z = xmin[3] + step*(iz-1)
        keep = false
        for iat in 1:length(solute)
          d2 = ( mysim.atoms[solute[iat]].coor[1] - x )^2 + 
               ( mysim.atoms[solute[iat]].coor[2] - y )^2 + 
               ( mysim.atoms[solute[iat]].coor[3] - z )^2  
          if d2 <= dmax2 && d2 > dmin2
            if keep == false
              keep = true
              igrid = igrid + 1
              grid[igrid] = MutableGMD3DGrid()
              grid[igrid].x = [ x, y, z ]
              grid[igrid].atom = solute[iat]
              grid[igrid].dmin = sqrt(d2)
            end
            if d2 <= grid[igrid].dmin
              grid[igrid].atom = solute[iat]
              grid[igrid].dmin = sqrt(d2)
            end
          end
        end
      end
    end
  end

  return grid

end

