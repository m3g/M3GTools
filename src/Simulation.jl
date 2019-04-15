include("./LogData.jl")

struct Simulation

  psf :: String
  dcd :: String

  natoms :: Int64
  atom :: Vector{Atom}

  nframes :: Int64
  dcdaxis :: Bool

  vmd :: String

  FortranDCD :: FortranFile

  logfile :: String
  log :: LogData

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

