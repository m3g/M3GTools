"""

In this example, we compute the time correlation function of the
rotations of the side chain of a trp residue

"""

#push!(LOAD_PATH,"../src")
#using Namd
include("../src/Namd.jl")

println(" Loading Plots... ")
using Plots

# Initialize simulation data and path to VMD executable
println(" Reading simulation data ... ")
Namd.init(psf="./structure.psf",
          dcd="./structure.dcd",
          vmd="vmd")

println(" Defining selections... ")

# Absorption vector, start and end
abs_start = Namd.select("index 1548 1549")
abs_end = Namd.select("index 1546")

# Absorption vector, start and end
emi_start = Namd.select("index 1548 1549")
emi_end = Namd.select("index 1545")

t, legendre, tcf = Namd.tcf(abs_start,abs_end,emi_start,emi_end;
                            r0=0.4,
                            theta=34.24,
                            scaletime=0.001)

println(" Plotting... ")
plot(t,legendre)
savefig("tcf.pdf")




