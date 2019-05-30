include("./LogData.jl")

struct Simulation

  psf :: String
  pdb :: String
  dcd :: String

  natoms :: Int64
  atoms :: Vector{Atom}

  nbonds :: Int64
  bonds :: Vector{Bond}

  nanbles :: Int64
  angles :: Vector{Angle}

  ndihedrals :: Int64
  dihedrals :: Vector{Dihedral}

  nimpropers :: Int64
  improper :: Vector{Improper}

  nframes :: Int64
  dcdaxis :: Bool

  vmd :: String

  FortranDCD :: FortranFile

  logfile :: String
  log :: LogData

  parfiles :: Vector{String}

end

function Base.show( io :: IO, simulation :: Simulation )
  println(" Namd simulation with ", simulation.natoms," atoms. ")
  if simulation.logfile != nothing
    println("   LOG file: ", simulation.logfile ) 
  end
  if simulation.psf != nothing
    println("   PSF file: ", simulation.psf ) 
  end
  if simulation.psf != nothing
    println("   DCD file: ", simulation.dcd ) 
  end
end

