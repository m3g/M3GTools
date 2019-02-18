
"""

In this example, we plot the temperature, read from the log file, as function
of simulation time

"""

push!(LOAD_PATH,"../")
using Namd

println(" Loading Plots... ")
using Plots
ENV["GKSwstype"]="nul" # This supresses the need of a display while plotting

# Initialize simulation data and path to VMD executable
println(" Reading simulation data ... ")
mysim = Namd.init(psf="./structure.psf",
                  log="./example.log")

println(" Plotting... ")
plot(mysim.log.time,mysim.log.temperature,label="example.log")
plot!(xlabel="time / ns")
plot!(ylabel="temperature / K")

savefig("temperature.pdf")


