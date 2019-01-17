
module Namd

  using FortranFiles

  # 
  # Scalars of the simulation
  # 

  local natoms = 0
  local nframes = 0
  local dcdaxis = true
  local psffile = "none"
  local dcdfile = "none"
  local mass :: Vector{}
  local charge :: Vector{}
  vmd_exec="vmd"

  #struct Atom
  #  index ::  Int64
  #  mass :: Float64
  #  charge :: Float64
  #  type :: FString{4}
  #  class :: FString{4}
  #  resname :: FString{4}
  #  resid :: Int64
  #  segname :: Fstring{4}
  #end 
  #atom(;index=0,mass=0.,charge=0.,type="X",
  #      class="X,resname="XXX",resid=0,segname="XXXX") =
  #     Atom(index,mass,charge,type,class,resname,resid,segname)
  
  # 
  # Initialize simulation data
  #

  function init(;psf=psfile,
                 dcd=dcdfile,
                 vmd=vmd_exec)

    global psffile = psf
    global dcdfile = FortranFile(dcd)
    global vmd_exec = vmd

    #
    # Read number of atoms and masses from PSF file
    #

    file = Base.open(psffile)
    start_atoms = false
    iatom = 0
    for line in eachline(file)
      data = split(line)
      if start_atoms 
        iatom = iatom + 1
        mass[iatom] = parse(Float32,data[8])
        charge[iatom] = parse(Float32,data[7])
        if iatom == natoms 
          break
        end
      end 
      if length(data) > 1 
        if data[2] == "!NATOM"
          global natoms = parse(Int64,data[1]) 
          global mass = Vector{Float32}(undef,natoms)
          global charge = Vector{Float32}(undef,natoms)
          start_atoms = true
          println(" Number of atoms: ", natoms)
        end 
      end 
    end 
    Base.close(file)
    println(" Masses were read to Namd.mass vector ")
    println(" Charges were read to Namd.charge vector ")

    #
    # Reads DCD file header, returns nframes (correctly, if set) and ntotat
    #

    IntVec = Vector{Int32}(undef,17)
    hdr, read_nframes, IntVec[1:8], ts, IntVec[9:17] = read(dcdfile, FString{4}, Int32, (Int32,8), Float64, (Int32,9))
    dummyi, title = read(dcdfile, Int32, FString{80})
    read_natoms = read(dcdfile,Int32) 

    if read_natoms != natoms
      error("ERROR: Number of atoms in PSF file differs from that in DCD file.")
    end
    global nframes = read_nframes
    println(" Number of frames in DCD file: ", nframes)

    # Check if dcd file contains axis information
 
    x = 0.
    try
      x = read(dcdfile, [ Float32 for i in 1:natoms ])
      dcdaxis = false
    catch err
      dcdaxis = true
    end
    if dcdaxis 
      println(" DCD contains periodic cell data. ")
    else
      println(" WARNING: DCD does NOT contain periodic cell data! ")
    end

    # Rewind the file and jump to first frame again
    rewind(dcdfile)
    for i in 1:3
      read(dcdfile)
    end

  end

  # 
  # Select atoms using vmd selection syntax, with vmd in background
  #

  function select(selection;update=false,
                            sides::Vector{Float64},
                            x::Vector{Float32},
                            y::Vector{Float32},
                            z::Vector{Float32})

    index_list = String
    readnext = false

    # If the selection has to be updated at every frame, we need
    # to write a temporary DCD file here to VMD

    if update
      writedcd(natoms,1,dcdaxis,sides,x,y,z;filename="Namdjl_DCDTEMP.dcd")
    end

    vmd_input = Base.open("./VMDINPUT_TMP.VMD","w")
    Base.write(vmd_input,"mol new \"$psffile\" \n")
    Base.write(vmd_input,"set sel [ atomselect top \"$selection\" ] \n")
    Base.write(vmd_input,"puts \"INDEXLIST\" \n")
    Base.write(vmd_input,"set indexes [ \$sel get index ] \n")
    Base.write(vmd_input,"puts \"ENDINDEXLIST\" \n")
    Base.write(vmd_input,"exit \n")
    Base.close(vmd_input)

    vmd_output = read(`$vmd_exec -dispdev text -e ./VMDINPUT_TMP.VMD`, String)

    for line in split(vmd_output,"\n")
      if readnext
        if line == "ENDINDEXLIST" 
          error("ERROR: Selection '$selection' does not contain any atom")
        end 
        index_list = line
        break
      end
      if line == "INDEXLIST" 
        readnext = true
      end
    end
    index_split = split(index_list)
    nsel = length(index_split)
    selection_indexes = Vector{Int64}(undef,nsel) 
    for i in 1:nsel
      selection_indexes[i] = parse(Int64,index_split[i]) + 1
    end

    run(`\rm -f ./VMDINPUT_TMP.VMD`)

    if ! update 
      println(" Selection '$selection' contains ",nsel," atoms ")
    end
    return selection_indexes

  end

  #
  # Fuctions to read the DCD file
  #

  closedcd() = close(dcdfile)
  rewinddcd() = rewind(dcdfile)

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

  function nextframe(;lastatom=0)

    if natoms == 0
      error("ERROR: Set natoms before reading frames with namd.header()")
    else
      lastatom = natoms
    end
    
    if dcdaxis 
      sides_read = read(dcdfile,(Float64,6))
    end
    x = read(dcdfile,(Float32,lastatom))
    y = read(dcdfile,(Float32,lastatom))
    z = read(dcdfile,(Float32,lastatom))
    
    sides = [ sides_read[1], sides_read[3], sides_read[6] ]

    return sides, x, y, z

  end

  #
  # Function that writes a dcd file
  #

  function writedcd(natoms,nframes,dcdaxis,sides,x,y,z;filename="Namdjl_DCDTEMP.dcd")

    dcdtemp = FortranFile(filename,"w";marker=RECMRK4B)

    IntVec = Vector{Int32}(undef,17)

    #hdr, read_nframes, IntVec[1:8], ts, IntVec[9:17] = 
    #  read(dcdfile, FString{4}, Int32, (Int32,8), Float64, (Int32,9))
    #dummyi, title = read(dcdfile, Int32, FString{80})
    #read_natoms = read(dcdfile,Int32) 

    hdr = "CORD"
    IntVec = [0 for i in 1:17]
    ts = 0.

#voltar: não funcionou
    write(dcdtemp,FString(4,hdr),Int32(nframes),IntVec[1:8],Float64(ts),IntVec[9:17])
    write(dcdtemp,Int32(1),FString(80,"DCDTEMP created by Namdjl"))
    write(dcdtemp,Int32(natoms))

    if nframes == 1 
      if dcdaxis
        write(dcdtemp,Float64(sides[1]),0.,Float64(sides[2]),0.,0.,Float64(sides[3]))
      end
      write(dcdtemp,x)
      write(dcdtemp,y)
      write(dcdtemp,z)
    #else
    #  for i in 1:nframes
    #    if dcdaxis
    #      write(dcdtemp,Float64(sides[i,1]),0.,Float64(sides[i,2]),0.,0.,Float64(sides[i,3]))
    #    end
    #    write(dcdtemp,(x[i,j],j=1,natoms))
    #    write(dcdtemp,(y[i,j],j=1,natoms))
    #    write(dcdtemp,(z[i,j],j=1,natoms))
    #  end
    end
    close(dcdtemp)

  end

  include("./cm.jl")
  include("./pbc.jl")
  include("./distance.jl")
  include("./elec.jl")
  include("./tcf.jl")
  include("./procrustes.jl")

end

