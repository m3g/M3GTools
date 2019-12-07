#
# Function that fills the grid value for each point of the read by reading the
# the file containing solute contributions 
#

function gmd3D(mysim :: Simulation, solute :: Vector{Int64}, gmd_solute_file :: String;
               dmax = 7.0, step = 0.5, dmin = 1.0)
            
  # Compute the grid points and get the closest atom to each grid point

  grid = gmd3D_setgrid(mysim,solute,dmax=dmax,step=step,dmin=dmin)

  # Now will read the solute_contribution file

  gmd_solute = readdlm(gmd_solute_file,comments=true,comment_char='#')

  # We need to get the atom indexes from the header of the file, to see to what column each atom
  # corresponds

  icolumn = Vector{Int64}(undef,length(solute))

  legacy = false
  file = open(gmd_solute_file,"r")
  for line in eachline(file)
    if line[1:1] != "#"
      break
    end
    if findfirst("mass:",line) == nothing
      continue
    end
    data = split(line)
    iat = 0
    try
      iat = parse(Int64,data[3])
    catch
      iat = parse(Int64,data[2])
      legacy = true
    end
    icol = parse(Int64,data[2]) + 2
    icolumn[iat] = icol
  end
  close(file)
  if legacy
    println(" Warning: Reading legacy GMD output (previous to 19.340) ")
    println("          if the atoms of the solute are not the first atoms, you may have problems.")
  end
  
  # Now, what we need to do is to, for each point in the grid, find the corresponding
  # atom in the gmd_solute table, and interpolate the gmd provided at the distance dmin

  ngrid = length(grid)

  for igrid in 1:ngrid           
    # Find values
    icol = icolumn[grid[igrid].atom]
    i = 1
    rup = gmd_solute[i,1]
    dup = gmd_solute[i,icol]
    rdown = 0.
    ddown = 0.
    while rup < grid[igrid].dmin
      rdown = rup
      ddown = dup
      i = i + 1
      rup = gmd_solute[i,1]
      dup = gmd_solute[i,icol]
    end
    # Interpolating
    grid[igrid].rho = ddown + ((dup-ddown)/(rup-rdown))*(grid[igrid].dmin-rdown)
  end

  # Now, returning the immutable array with the final grid

  gmd3D = Vector{GMD3DGrid}(undef,ngrid)
  for i in 1:ngrid
    gmd3D[i] = GMD3DGrid(grid[i].x,grid[i].atom,grid[i].dmin,grid[i].rho)
  end
  
  return filter( x -> x.rho > 0. , gmd3D )

end
