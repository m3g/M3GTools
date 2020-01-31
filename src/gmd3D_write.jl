#
# Function that writes the gmd3D to a pdb file
#

# If called already with the grid of the gmd3D computed

function gmd3D_write(mysim :: Simulation, grid :: Vector{GMD3DGrid}, output :: String; scale=nothing)

  n = length(grid)

  if scale == nothing
    minrho = grid[1].rho
    maxrho = grid[1].rho
    for i in 2:n
      minrho = min(minrho,grid[i].rho)
      maxrho = max(maxrho,grid[i].rho)
    end
  else
    minrho = scale[1]
    maxrho = scale[2]
  end

  file = open(output,"w")
  println(file,"REMARK  3D representation of MDDF density ")
  println(file,"REMARK  ")
  println(file,"REMARK  Occupancy column contains distance to solute. ")
  println(file,"REMARK  B-factor column contains the scaled density, such that 99.9 is the maximum ")
  println(file,"REMARK  density obsered and 0. is the minimum density observed. The actual limits are: ")
  println(file,"REMARK  Minimum density: $minrho")
  println(file,"REMARK  Maximum density: $maxrho")
  for i in 1:n
    iat = grid[i].atom
    name = mysim.atoms[iat].name
    resname = mysim.atoms[iat].resname
    chain = mysim.atoms[iat].chain
    resnum = mysim.atoms[iat].resid
    x = grid[i].x[1]
    y = grid[i].x[2]
    z = grid[i].x[3]
    b = 99.9*(grid[i].rho - minrho)/(maxrho-minrho)
    occup = grid[i].dmin
    model = 0
    atom = PDBTools.Atom(i,name,resname,chain,resnum,x,y,z,b,occup,model)
    println(file,PDBTools.write_atom(atom))
  end
  close(file)

end

# If called directly with the solute selection and gmd_solute file name

function gmd3D_write(mysim :: Simulation, 
                     solute :: Vector{Int64}, 
                     gmd_solute_file :: String, 
                     output :: String;
                     dmax = 7.0, dmin = 1.0, step = 0.5, scale=nothing)

  grid = gmd3D(mysim, solute, gmd_solute_file, 
               dmax = 7.0, dmin = 1.0, step = 0.5)

  gmd3D_write(mysim, grid, output, scale=scale)

end










