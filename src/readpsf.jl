#
# Read PSF file
#

function readpsf(psf)

  #
  # Read number of atoms and masses from PSF file
  #

  local natoms, atom, dcdaxis

  file = Base.open(psf)
  start_atoms = false
  iatom = 0
  residue = 0
  for line in eachline(file)
    data = split(line)
    if start_atoms 
      iatom = iatom + 1
      segname = data[2]
      resid = parse(Int64,data[3])
      if iatom == 1 
        residue = 1
      else
        if atom[iatom-1].resid != resid
          residue = residue + 1
        end
      end
      resname = data[4]
      name = data[5]
      type = data[6]
      charge = parse(Float32,data[7])
      mass = parse(Float32,data[8])
      atom[iatom] = Atom(iatom,residue,resid,name,resname,segname,type,charge,mass)
      if iatom == natoms 
        break
      end
    end 
    if length(data) > 1 
      if data[2] == "!NATOM"
        natoms = parse(Int64,data[1]) 
        atom = Vector{Atom}(undef,natoms)
        start_atoms = true
      end 
    end 
  end 
  Base.close(file)

  return atom

end

