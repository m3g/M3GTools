
module Namd

  using FortranFiles

  # 
  # Structures
  # 

  struct LogData

   timestep :: Int32
   dcdfreq :: Int32
   dcdfirststep :: Int32
   
   time :: Vector{Float64}

   ts :: Vector{Int64}
   bond :: Vector{Float64}
   angle :: Vector{Float64}
   dihed :: Vector{Float64}
   imprp :: Vector{Float64}
   elect :: Vector{Float64}
   vdw :: Vector{Float64}
   boundary :: Vector{Float64}
   misc :: Vector{Float64}
   kinetic :: Vector{Float64}
   total :: Vector{Float64}
   temperature :: Vector{Float64}
   potential :: Vector{Float64}
   total3 :: Vector{Float64}
   tempavg :: Vector{Float64}
   pressure :: Vector{Float64}
   gpressure :: Vector{Float64}
   volume :: Vector{Float64}
   pressavg :: Vector{Float64}
   gpressavg :: Vector{Float64}

  end
  LogData() = LogData(0, 0, 0, [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], 
                      [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.])

  struct Simulation

    psf :: String
    dcd :: String

    natoms :: Int64
    nframes :: Int64
    dcdaxis :: Bool

    mass :: Vector{Float64}
    charge :: Vector{Float64}

    vmd :: String

    FortranDCD :: FortranFile
    log :: LogData

  end

  function init(;psf="none",
                 dcd="none",
                 log="none",
                 vmd="vmd")

    #
    # Read number of atoms and masses from PSF file
    #
  
    local natoms, mass, charge, dcdaxis

    if psf == "none"
      error(" At least a PSF must be provided with psf=filename.psf ")
    end

    file = Base.open(psf)
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
          natoms = parse(Int64,data[1]) 
          mass = Vector{Float32}(undef,natoms)
          charge = Vector{Float32}(undef,natoms)
          start_atoms = true
          println(" Number of atoms: ", natoms)
        end 
      end 
    end 
    Base.close(file)
    println(" Masses were read to mass vector ")
    println(" Charges were read to charge vector ")

    #
    # Reads DCD file header, returns nframes (correctly, if set) and ntotat
    #

    if dcd != "none"
      FortranDCD = FortranFile(dcd)
      nframes, read_natoms, dcdaxis = dcdheader(FortranDCD) 
      if read_natoms != natoms
        error(" Number of atoms in DCD file is different from that of PSF file.")
      end
    else
      nframes = 0
      dcdaxis = false
      FortranDCD = FortranFile(psf)
      println(" Warning: no DCD file provided, so many functions won't work. ")
    end

    #
    # Read data from log file, if provided
    #

    if log != "none" 

      println(" Reading data from LOG file: ", strip(log))

      timestep = 0
      dcdfreq = 0
      dcdfirststep = 0
      nsteps = 0
      for line in eachline(log)
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
      for line in eachline(log)
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
      logdata = LogData()
    end

    simulation = Simulation(psf,dcd,
                            natoms,nframes,dcdaxis,
                            mass,charge,
                            vmd,
                            FortranDCD,
                            logdata)

    return simulation

  end

  include("./select.jl")
  include("./dcdio.jl")

  include("./cm.jl")
  include("./pbc.jl")
  include("./distance.jl")
  include("./elec.jl")
  include("./tcf.jl")
  include("./procrustes.jl")

end

