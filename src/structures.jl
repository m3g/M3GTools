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

struct Atom

  index :: Int64
  residue :: Int64
  resid :: Int64
  
  name :: String
  resname :: String
  segname :: String
  type :: String

  charge :: Float64
  mass :: Float64
  
  backbone :: Bool

end

struct Bond

   index1 :: Int64
   index2 :: Int64

end

struct Angle

   index1 :: Int64
   index2 :: Int64
   index3 :: Int64

end

struct Dihedral

   index1 :: Int64
   index2 :: Int64
   index3 :: Int64
   index4 :: Int64

end

struct Improper

   index1 :: Int64
   index2 :: Int64
   index3 :: Int64
   index4 :: Int64

end

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
  if simulation.logfile != "none"
    println("   LOG file: ", simulation.logfile ) 
  end
  if simulation.psf != "none"
    println("   PSF file: ", simulation.psf ) 
  end
  if simulation.psf != "none"
    println("   DCD file: ", simulation.dcd ) 
  end
end

function Base.show( io :: IO, log :: LogData )
  ndata = length(log.time)
  tempavg = 0.
  pavg = 0.
  eavg = 0.
  for i in 1:ndata
    tempavg = tempavg + log.temperature[i]
    pavg = pavg + log.pressure[i]
    eavg = eavg + log.total[i]
  end
  tempavg = tempavg / ndata
  pavg = pavg / ndata
  eavg = eavg / ndata
  println(" Simulation log data: ")
  println("    Number of steps printed: ", ndata)
  println("    Average temperature: ", tempavg)
  println("    Average pressure: ", pavg)
  println("    Average total energy: ", eavg)
end

