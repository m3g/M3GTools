#
# Subroutine readpsf: Read all that from the psf file and allocates
#                     some arrays of the force_field module.
#
# L. Martinez, Nov 11, 2013
# Institute of Chemistry - State University of Campinas - Brazil
#
# Translated to Julia: April 14, 2019
#

function readpsf(psf)

  file = open(psf,"r")
  lines = readlines(file)
  close(file)
  
  # Get the number of atoms, bonds, etc, from the PSF file headers

  local natoms = 0
  local nbonds = 0 
  local nangles = 0 
  local ndihed = 0 
  local nimpr = 0 
  local ncrterm = 0

  for line in lines

    data = split(line)
    if length(data) < 2
      continue
    end

    if data[2] == "!NATOM" 
      natoms = parse(Int64,data[1])
    end
    if data[2] == "!NBOND:"
      nbonds = parse(Int64,data[1])
    end
    if data[2] == "!NTHETA:"
      nangles = parse(Int64,data[1])
    end
    if data[2] == "!NPHI:"
      ndihed = parse(Int64,data[1])
    end
    if data[2] == "!NIMPHI:"
      nimpr = parse(Int64,data[1])
    end
    if data[2] == "!NCRTERM:"
      ncrterm = parse(Int64,data[1])
    end
      
  end

  if natoms == 0 
    error(" ERROR: Flag !NATOM not found in PSF file. ")
  end
  
  # Allocate all arrays

  atoms = Vector{Atom}(undef,natoms)
  bonds = Vector{Bond}(undef,nbonds)
  angles = Vector{Angle}(undef,nangles)
  dihedrals = Vector{Dihedral}(undef,ndihed)
  impropers = Vector{Improper}(undef,nimpr)
  
  # Read PSF file data

  local residue

  local read_atoms = false
  local read_bonds = false
  local read_angles = false
  local read_dihed = false
  local read_impr = false

  local iatom = 0
  local ibond = 0
  local iangle = 0
  local idihed = 0
  local iimpr = 0

  for line in lines
  
    data = split(line)
    if length(data) < 2 
      continue
    end

    #
    # Read atoms
    #

    if data[2] ==  "!NATOM" 
      read_atoms = true
      continue
    end
    if read_atoms

      iatom = iatom + 1
      data = split(line)

      index = parse(Int64,data[1])
      if iatom == 1
        residue = 1
      else
        if parse(Int64,data[3]) != atoms[iatom-1].residue
          residue = residue + 1
        end
      end
      resid = parse(Int64,data[3])
      name = data[5]
      resname = data[4]
      segname = data[2]
      type = data[6]
      charge = parse(Float64,data[7])
      mass = parse(Float64,data[8])

      backbone = false
      if ( name in ["HN","N","HA","CA","C","O"] ) ||
         ( resname == "GLY" && name in ["HA1","HA2"] )
        backbone = true
      end

      atoms[iatom] = Atom(index,residue,resid,name,resname,segname,type,charge,mass,backbone)

      if iatom == natoms 
        read_atoms = false
        continue
      end
    end

    #
    # Read bonds
    #

    if data[2] == "!NBOND:"
      read_bonds = true
      continue
    end
    if read_bonds
      data = split(line)
      ind = @. parse(Int64,data)
      for i in 1:2:length(data)-1
        ibond = ibond + 1
        bonds[ibond] = Bond(ind[i],ind[i+1])
      end
      if ibond == nbonds
        read_bonds = false
        continue
      end
    end

    #
    # Read angles
    #

    if data[2] == "!NTHETA:"
      read_angles = true
      continue
    end
    if read_angles
      data = split(line)
      ind = @. parse(Int64,data)
      for i in 1:3:length(data)-1
        iangle = iangle + 1
        angles[iangle] = Angle(ind[i],ind[i+1],ind[i+2])
      end
      if iangle == nangles
        read_angles = false
        continue
      end
    end

    #
    # Read dihedrals
    #

    if data[2] == "!NPHI:"
      read_dihed = true
      continue
    end
    if read_dihed
      data = split(line)
      ind = @. parse(Int64,data)
      for i in 1:4:length(data)-1
        idihed = idihed + 1
        dihedrals[idihed] = Dihedral(ind[i],ind[i+1],ind[i+2],ind[i+3])
      end
      if idihed == ndihed
        read_dihed = false
        continue
      end
    end

    #
    # Read impropers
    #

    if data[2] == "!NIMPHI:"
      read_impr = true
      continue
    end
    if read_impr
      data = split(line)
      ind = @. parse(Int64,data)
      for i in 1:4:length(data)-1
        iimpr = iimpr + 1
        impropers[iimpr] = Improper(ind[i],ind[i+1],ind[i+2],ind[i+3])
      end
      if iimpr == nimpr
        read_impr = false
        continue
      end
    end

  end

  return atoms, bonds, angles, dihedrals, impropers
end

import Base.show
Base.show( io :: IO, atoms :: Vector{Atom} ) = print("$(length(atoms)) atoms")
Base.show( io :: IO, bonds :: Vector{Bond} ) =  print("$(length(bonds)) angles")
Base.show( io :: IO, angles :: Vector{Angle} ) =  print("$(length(angles)) bonds")
Base.show( io :: IO, dihedrals :: Vector{Dihedral} ) =  print("$(length(dihedrals)) dihedrals")
Base.show( io :: IO, impropers :: Vector{Improper} ) =  print("$(length(impropers)) improper dihedrals")

