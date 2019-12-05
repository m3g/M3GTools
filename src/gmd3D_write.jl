#
# Function that writes the gmd3D to a pdb file
#

# If called already with the gmd3D computed

function gmd3D_write(mysim :: Simulation, grid :: Vector{GMD3DGrid}, output :: String)

  n = length(gmd3D)
  minrho = gmd3D[1].rho
  maxrho = gmd3D[1].rho
  for i in 2:n
    minrho = min(minrho,gmd3D[i].rho)
    maxrho = max(maxrho,gmd3D[i].rho)
  end

  file = open(output,"w")
  for i in 1:n
    iat = gmd3D[i].atom
    name = mysim.atoms[iat].name
    resname = mysim.atoms[iat].resname
    chain = mysim.atoms[iat].chain
    resnum = mysim.atoms[iat].resid
    x = gmd3D[i].x[1]
    y = gmd3D[i].x[2]
    z = gmd3D[i].x[3]
    b = 99.9*(gmd3D[i].rho - minrho)/(maxrho-minrho)
    occup = gmd3D[i].dmin
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
                     dmax = 7.0, dmin = 1.0, step = 0.5)

  grid = gmd3D(mysim, solute, gmd_solute_file, 
               dmax = 7.0, dmin = 1.0, step = 0.5)

  gmd3D_write(mysim, grid, output)

end










