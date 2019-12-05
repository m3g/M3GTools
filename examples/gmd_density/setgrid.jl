#
# This function builds the grid of the 3D density function and fills an array of
# mutable structures of type MutableDensity, containing the position of the atoms of 
# grid, the closest atom to that position, and distance. 
#

function setgrid(mysim,solute;dmax=7.0,step=0.5,dmin=1.0)

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

  # Initializing density structure

  density = Vector{MutableDensity}(undef,ngrid)
 
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
              density[igrid] = MutableDensity()
              density[igrid].x = [ x, y, z ]
              density[igrid].atom = solute[iat]
              density[igrid].dmin = sqrt(d2)
            end
            if d2 <= density[igrid].dmin
              density[igrid].atom = solute[iat]
              density[igrid].dmin = sqrt(d2)
            end
          end
        end
      end
    end
  end

  return density

end

#
# Function that fills the density value for each point of the read by reading the
# the file containing solute contributions 
#

function density3D(mysim,solute,gmd_solute_file;
                   dmax = 7.0, step = 0.5, dmin = 1.0)
            
  # Compute the grid points and get the closest atom to each grid point

  density = setgrid(mysim,solute,dmax=dmax,step=step,dmin=dmin)

  # Now will read the solute_contribution file

  gmd_solute = readdlm(gmd_solute_file,comments=true,comment_char='#')

  # We need to get the atom indexes from the header of the file, to see to what column each atom
  # corresponds

  icolumn = Vector{Int64}(undef,length(solute))

  file = open(gmd_solute_file,"r")
  for line in eachline(file)
    if line[1:1] != "#"
      break
    end
    if findfirst("mass:",line) == nothing
      continue
    end
    data = split(line)
    iat = parse(Int64,data[3])
    icol = parse(Int64,data[2]) + 2
    icolumn[iat] = icol
  end
  close(file)
  
  # Now, what we need to do is to, for each point in the grid, find the corresponding
  # atom in the gmd_solute table, and interpolate the gmd provided at the distance dmin

  ngrid = length(density)

  for igrid in 1:ngrid           
    # Find values
    icol = icolumn[density[igrid].atom]
    i = 1
    rup = gmd_solute[i,1]
    dup = gmd_solute[i,icol]
    rdown = 0.
    ddown = 0.
    while rup < density[igrid].dmin
      rdown = rup
      ddown = dup
      i = i + 1
      rup = gmd_solute[i,1]
      dup = gmd_solute[i,icol]
    end
    # Interpolating
    density[igrid].rho = ddown + ((dup-ddown)/(rup-rdown))*(density[igrid].dmin-rdown)
  end

  # Now, returning the immutable array with the final grid

  density3D = Vector{Density}(undef,ngrid)
  for i in 1:ngrid
    density3D[i] = Density(density[i].x,density[i].atom,density[i].dmin,density[i].rho)
  end
  
  return filter( x -> x.rho > 0. , density3D )

end

using PDBTools # From https://github.com/mcubeg/PDBTools

function write_density3D(mysim,density3D,file)

  open(file,"w")
  n = length(density3D)
  for i in 1:n
    iat = density3D[i].atom
    name = mysim.atoms[iat].name
    resname = mysim.atoms[iat].resname
    chain = mysim.atoms[iat].chain
    resnum = mysim.atoms[iat].resid
    x = density3D[i].x[1]
    y = density3D[i].x[2]
    z = density3D[i].x[3]
    b = density3D[i].rho
    occup = density3D[i].dmin
    model = 0
    atom = PDBTools.Atom(i,name,resname,chain,resnum,x,y,z,b,occup,model)
    println(file,PDBTools.write_atom(atom))
  end
  close(file)

end










