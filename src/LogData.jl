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

