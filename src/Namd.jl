
module Namd

  using FortranFiles

  # 
  # Scalars of the simulation
  # 

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

  end

  function init(;psf="none",
                 dcd="none",
                 vmd="vmd")

    #
    # Read number of atoms and masses from PSF file
    #
  
    local natoms, mass, charge, dcdaxis

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

    FortranDCD = FortranFile(dcd)
    nframes, read_natoms, dcdaxis = dcdheader(FortranDCD) 
    if read_natoms != natoms
      error(" Number of atoms in DCD file is different from that of PSF file.")
    end

    simulation = Simulation(psf,dcd,
                            natoms,nframes,dcdaxis,
                            mass,charge,
                            vmd,
                            FortranDCD)

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

