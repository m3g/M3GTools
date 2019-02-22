#
# Functions to read the DCD file
#

using FortranFiles

# 
# Functions to close and rewind a dcd file
#

closedcd(simulation :: Simulation) = close(simulation.FortranDCD)

rewinddcd(simulation :: Simulation) = rewind(simulation.FortranDCD)

#
# Function to read the header given the simulation data
#

function dcdheader(simulation :: Simulation) 

  read_nframes, read_natoms, dcdaxis = dcdheader(simulation.FortranDCD)

  if read_natoms != simulation.natoms
    error("ERROR: Number of atoms in PSF file differs from that in DCD file.")
  end
  nframes = read_nframes
  println(" Number of frames in DCD file: ", nframes)

  if dcdaxis
    println(" DCD contains periodic cell data. ")
  else
    println(" WARNING: DCD does NOT contain periodic cell data! ")
  end

  return read_nframes, read_natoms, dcdaxis

end 

#
# Function to read the header given only the file stream
#

function dcdheader(FortranDCD :: FortranFile)

    rewind(FortranDCD)

    # Read header
    IntVec = Vector{Int32}(undef,17)
    hdr, read_nframes, IntVec[1:8], ts, IntVec[9:17] = read(FortranDCD, FString{4}, Int32, (Int32,8), Float64, (Int32,9))
    dummyi, title = read(FortranDCD, Int32, FString{80})
    read_natoms = read(FortranDCD,Int32)

    # Check if dcd file contains axis information
    dcdaxis = false
    x = 0.
    try
      x = read(FortranDCD, [ Float32 for i in 1:read_natoms ])
    catch err
      dcdaxis = true
    end

    firstframe(FortranDCD)

    return read_nframes, read_natoms, dcdaxis

end

# Function that goes to the file point to read the first frame

firstframe(simulation :: Simulation) = firstframe(simulation.FortranDCD)
function firstframe(FortranDCD :: FortranFile)

    # rewind
    rewind(FortranDCD)

    # skip header
    read(FortranDCD)
    read(FortranDCD)
    read(FortranDCD)

end

#
# Function to read next frame
#

function nextframe(simulation :: Simulation;lastatom=0)  

  if simulation.natoms == 0
    error(" ERROR: Set natoms before reading frames with Namd.init ")
  else
    lastatom = simulation.natoms
  end
  
  sides, x, y, z = nextframe(simulation.FortranDCD, simulation.dcdaxis, lastatom)
  return sides, x, y, z

end

function nextframe(FortranDCD :: FortranFile, dcdaxis, lastatom)

  if dcdaxis 
    sides_read = read(FortranDCD,(Float64,6))
    sides = [ sides_read[1], sides_read[3], sides_read[6] ]
  else
    sides = zeros(3)
  end 
  x = read(FortranDCD,(Float32,lastatom))
  y = read(FortranDCD,(Float32,lastatom))
  z = read(FortranDCD,(Float32,lastatom))
  
  return sides, x, y, z

end

#
# Function to read the number of frames if the header does not inform it
#

function getnframes(FortranDCD :: FortranFile, dcdaxis :: Bool )
  firstframe(FortranDCD)
  nframes = 0
println(dcdaxis)
  while true
    try 
      if dcdaxis
        x = read(FortranDCD,Float64)
        println(x)
      end
      x = read(FortranDCD,Float32)
        println(x)
      x = read(FortranDCD,Float32)
        println(x)
      x = read(FortranDCD,Float32)
        println(x)
      nframes = nframes + 1
    catch
      firstframe(FortranDCD)
      return nframes
    end
  end
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

