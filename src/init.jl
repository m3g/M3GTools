#
# Initialize simulation data
#

function init(;psf=nothing,
               dcd=nothing,
               log=nothing,
               pdb=nothing,
               parfiles=nothing,
               vmd="vmd")

  #
  # Read PSF file data
  #

  if psf == nothing
    error(" At least a PSF must be provided with psf=filename.psf ")
  end
  atoms, bonds, angles, dihedrals, impropers = readpsf(psf)
  natoms = length(atoms)
  nbonds = length(bonds)
  nangles = length(angles)
  ndihedrals = length(dihedrals)
  nimpropers = length(impropers)

  # Reading coordinates from PDB file, if provided

  if pdb != nothing
    getpdbcoords!(pdb,atoms)
  else
    pdb = "none"
  end

  # Reading parameter files, if provided

  if parfiles != nothing
    readprm!(parfiles,atoms,bonds,angles,dihedrals,impropers)
  else
    parfiles = ["none"]
  end

  #
  # Reads DCD file header, returns nframes (correctly, if set) and ntotat
  #

  if dcd != nothing
    FortranDCD = FortranFile(dcd)
    nframes, read_natoms, dcdaxis = dcdheader(FortranDCD) 
    if read_natoms != natoms
      error(" Number of atoms in DCD file is different from that of PSF file.")
    end
    if nframes < 1
      println(" WARNING: The DCD header does not inform the number of frames. Will read the DCD to fetch it. ")
      nframes = getnframes(FortranDCD, dcdaxis)
    end
  else
    dcd = "none"
    nframes = 0
    dcdaxis = false
    FortranDCD = FortranFile(psf)
    println(" Warning: no DCD file provided, so many functions won't work. ")
  end

  #
  # Read data from log file, if provided
  #

  if log != nothing

    logfile = strip(log)
    println(" Reading data from LOG file: ", logfile,"\n" )

    timestep = 0
    dcdfreq = 0
    dcdfirststep = 0
    nsteps = 0
    for line in eachline(logfile)
      if line[1:min(length(line),7)] == "ENERGY:"
        nsteps = nsteps + 1
      end
      try 
        if line[1:14] == "Info: TIMESTEP"
          data = split(line)
          timestep = parse(Int64,data[3])
        end
      catch msg
      end
      try 
        if line[1:19] == "Info: DCD FREQUENCY"
          data = split(line)
          dcdfreq = parse(Int64,data[4])
        end
      catch msg
      end
      try 
        if line[1:19] == "Info: DCD FIRST STEP"
          data = split(line)
          dcdfirststep = parse(Int64,data[5])
        end
      catch msg
      end
    end
    if timestep == 0 
      error("Could not find TIMESTEP in log file.")
    end

    time = Vector{Float64}(undef,nsteps)
    ts = Vector{Int64}(undef,nsteps)
    bond = Vector{Float64}(undef,nsteps)
    angle = Vector{Float64}(undef,nsteps)
    dihed = Vector{Float64}(undef,nsteps)
    imprp = Vector{Float64}(undef,nsteps) 
    elect = Vector{Float64}(undef,nsteps) 
    vdw = Vector{Float64}(undef,nsteps) 
    boundary = Vector{Float64}(undef,nsteps) 
    misc = Vector{Float64}(undef,nsteps) 
    kinetic = Vector{Float64}(undef,nsteps) 
    total = Vector{Float64}(undef,nsteps) 
    temperature = Vector{Float64}(undef,nsteps) 
    potential = Vector{Float64}(undef,nsteps) 
    total3 = Vector{Float64}(undef,nsteps) 
    tempavg = Vector{Float64}(undef,nsteps) 
    pressure = Vector{Float64}(undef,nsteps) 
    gpressure = Vector{Float64}(undef,nsteps) 
    volume = Vector{Float64}(undef,nsteps) 
    pressavg = Vector{Float64}(undef,nsteps) 
    gpressavg = Vector{Float64}(undef,nsteps) 

    nsteps = 0
    for line in eachline(logfile)
      if line[1:min(length(line),7)] == "ENERGY:"
        nsteps = nsteps + 1
        data = split(line)
        ts[nsteps] = parse(Int64,data[2])
        time[nsteps] = ts[nsteps] / ( 1.e6 * timestep ) # in ns
        bond[nsteps] = parse(Float64,data[3])
        angle[nsteps] = parse(Float64,data[4]) 
        dihed[nsteps] = parse(Float64,data[5]) 
        imprp[nsteps] = parse(Float64,data[6]) 
        elect[nsteps] = parse(Float64,data[7]) 
        vdw[nsteps] = parse(Float64,data[8]) 
        boundary[nsteps] = parse(Float64,data[9]) 
        misc[nsteps] = parse(Float64,data[10]) 
        kinetic[nsteps] = parse(Float64,data[11]) 
        total[nsteps] = parse(Float64,data[12]) 
        temperature[nsteps] = parse(Float64,data[13]) 
        potential[nsteps] = parse(Float64,data[14]) 
        total3[nsteps] = parse(Float64,data[15]) 
        tempavg[nsteps] = parse(Float64,data[16]) 
        pressure[nsteps] = parse(Float64,data[17]) 
        gpressure[nsteps] = parse(Float64,data[18]) 
        volume[nsteps] = parse(Float64,data[19]) 
        pressavg[nsteps] = parse(Float64,data[20]) 
        gpressavg[nsteps] = parse(Float64,data[21]) 
      end
    end  
    logdata = LogData(timestep,dcdfreq,dcdfirststep,time,ts,bond,angle,dihed,imprp,
                      elect,vdw,boundary,misc,kinetic,total,temperature,potential,total3,tempavg,
                      pressure,gpressure,volume,pressavg,gpressavg)
  else
    logfile = "none"
    logdata = LogData()
  end

  simulation = Simulation(psf,pdb,dcd,
                          natoms,atoms,
                          nbonds,bonds,
                          nangles,angles,
                          ndihedrals,dihedrals,
                          nimpropers,impropers,
                          nframes,dcdaxis,
                          vmd,
                          FortranDCD,
                          logfile,logdata,parfiles)

  return simulation
end

