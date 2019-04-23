#
# Function to read atomic coordinates from a PDB file.
#
function getpdbcoords!(pdbfile,atoms)
  file = open(pdbfile,"r")
  lines = readlines(file)
  local iatom = 0
  for line in lines
    data = split(line)
    if length(data) < 1
      continue
    end
    if data[1] == "ATOM" || data[1] == "HETATM"
      iatom = iatom + 1
      coor = [ parse(Float64,line[31:38]), parse(Float64,line[39:46]),  parse(Float64,line[47:54]) ]
      a = atoms[iatom]
      atoms[iatom] = Atom(a.index,a.residue,a.resid,a.name,a.resname,a.segname,a.type,a.element,
                          a.charge,a.mass,a.eps,a.rmin,a.eps14,a.rmin14,a.backbone,
                          coor)
    end
  end
end
