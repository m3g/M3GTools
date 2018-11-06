
module namd

  using FortranFiles

  # 
  # Scalars of the simulation
  # 

  natoms = 0
  nframes = 0
  dcdaxis = true
  
  # 
  # Set name of PSF file
  #

  psffile="none"
  function psf(user_psf_file) 
    global psffile = user_psf_file
  end
  function checkpsf(psffile)
   if psffile == "none"
     error("ERROR: You must define the psf file first, with: namd.psf(\"file.psf\")")
   end
  end

  #
  # Get number of atoms from psf file 
  #
  
  function psf_natoms(;psffile=psffile)

    checkpsf(psffile)

    file = Base.open(psffile)
    natoms = 0
    for line in eachline(file)
      data = split(line)
      if length(data) > 1 
        if data[2] == "!NATOM"
          natoms = parse(Int64,data[1]) 
          break
        end 
      end 
    end 
    Base.close(file)

    return natoms

  end

  # 
  # Read masses from PSF file 
  #

  function mass(;psffile=psffile)

    checkpsf(psffile)
    mass = Vector{}
    natoms :: Int64 = 0

    file = Base.open(psffile)
    start_atoms = false
    iatom = 0
    for line in eachline(file)
      data = split(line)
      if start_atoms 
        iatom = iatom + 1
        mass[iatom] = parse(Float32,data[8])
        if iatom == natoms 
          break
        end
      end 
      if length(data) > 1 
        if data[2] == "!NATOM"
          natoms = parse(Int64,data[1]) 
          mass = Vector{Float32}(undef,natoms)
          start_atoms = true
        end 
      end 
    end 
    Base.close(file)

    return mass

  end
  
  #
  # Set vmd executable
  #

  vmd_exec = "vmd"
  function vmd(user_vmd_exec)
    global vmd_exec = user_vmd_exec
  end

  # 
  # Select atoms using vmd selection syntax, with vmd in background
  #

  function select(selection;psffile=psffile,vmd_exec=vmd_exec)

    checkpsf(psffile)

    index_list = String
    readnext = false

    vmd_input = Base.open("./VMDINPUT_TMP.VMD","w")
    Base.write(vmd_input,"mol new \"$psffile\" \n")
    Base.write(vmd_input,"set sel [ atomselect top \"$selection\" ] \n")
    Base.write(vmd_input,"puts \"INDEXLIST\" \n")
    Base.write(vmd_input,"set indexes [ \$sel get index ] \n")
    Base.write(vmd_input,"exit \n")
    Base.close(vmd_input)

    vmd_output = read(`$vmd_exec -dispdev text -e ./VMDINPUT_TMP.VMD`, String)

    for line in split(vmd_output,"\n")
      if readnext
        index_list = line
        break
      end
      if line == "INDEXLIST" 
        readnext = true
      end
    end
    index_split = split(index_list)
    nsel = length(index_split)
    selection = Vector{Int64}(undef,nsel) 
    for i in 1:nsel
      selection[i] = parse(Int64,index_split[i]) + 1
    end

    run(`\rm -f ./VMDINPUT_TMP.VMD`)

    return nsel, selection

  end

  #
  # Read DCD files
  #

  closedcd() = close(dcdfile)
  rewinddcd() = rewind(dcdfile)

  #
  # Set name of dcd file 
  # 

  function checkdcd(dcdfile)
   if dcdfile == "none"
     error("ERROR: You must define the DCD file first, with: namd.dcd(\"file.dcd\")")
   end
  end

  #
  # Reads DCD file header, returns nframes (correctly, if set) and ntotat
  #

  dcdfile="none"
  function dcd(user_dcdfile)

    global dcdfile = FortranFile(user_dcdfile)
   
    Dint1 = Vector{Int32}(undef, 8)
    Dint2 = Vector{Int32}(undef, 9)
    dummyc, read_nframes, Dint1, dummyr, Dint2 = read(dcdfile, FString{4}, Int32, (Int32,8), Float32, (Int32,9))
    dummyi, dummyr = read(dcdfile, Int32, Float32)
    read_natoms = read(dcdfile,Int32)

    global natoms = read_natoms
    global nframes = read_nframes

    # Check if dcd file contains axis information
 
    x = 0.
    try
      x = read(dcdfile, [ Float32 for i in 1:natoms ])
      dcdaxis = false
    catch err
      dcdaxis = true
    end

    rewind(dcdfile)
    for i in 1:3
      read(dcdfile)
    end

  end

  #function nframes(dcdfile, dcdaxis, lastframe) 
  #  checkdcd(dcdfile)
  #  rewind(dcdfile)
  #  nframes, ntotat = header(dcdfile)
  #  nframes = 0
  #  
  #  if dcdaxis 
  #    read(file)
  #  end
  #  
  #  return nframes
  #end

  function nextframe(;dcdfile=dcdfile)

    if natoms == 0
      error("ERROR: Set natoms before reading frames with namd.header()")
    end
    
    sides_read = read(dcdfile,(Float64,6))
    x = read(dcdfile,(Float32,natoms))
    y = read(dcdfile,(Float32,natoms))
    z = read(dcdfile,(Float32,natoms))
    
    sides = [ sides_read[1], sides_read[3], sides_read[6] ]

    return sides, x, y, z

  end
  
end

function cm(selection,mass,x,y,z)

  cm = [ 0., 0., 0. ]
  totmass = 0.
  for i in 1:length(selection)
    cm[1] = cm[1] + mass[selection[i]]*x[selection[i]] 
    cm[2] = cm[2] + mass[selection[i]]*y[selection[i]] 
    cm[3] = cm[3] + mass[selection[i]]*z[selection[i]] 
    totmass = totmass + mass[selection[i]]
  end
  cm[1] = cm[1] / totmass
  cm[2] = cm[2] / totmass
  cm[3] = cm[3] / totmass

  return cm

end

