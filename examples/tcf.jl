"""

In this example, we compute the time correlation function of the
rotations of the side chain of a trp residue

"""

push!(LOAD_PATH,"../src")
using Namd

println(" Loading Plots... ")
using Plots
ENV["GKSwstype"]="nul" # This supresses the need of a display while plotting

# Initialize simulation data and path to VMD executable
println(" Reading simulation data ... ")
mysim = Namd.init(psf="./structure.psf",
                  dcd="./structure.dcd",
                  vmd="vmd")

println(" Defining selections... ")

# Absorption vector, start and end (the center of mass of selected atoms is considered)

abs_start = Namd.select(mysim,"index 1548 1549")
abs_end = Namd.select(mysim,"index 1546")

# Emission vector, start and end
#emi_start = Namd.select("index 1548 1549")
emi_start = Namd.select(mysim,"index 1548 1549")
emi_end = Namd.select(mysim,"index 1545")

# Optional: align coordinates to first frame, according to the alignment of
# the selection defined here. In this case, we remove the rotation of the protein.

CAs = Namd.select(mysim,"protein and name CA")

t, legendre, tcf = Namd.tcf(mysim,
                            abs_start,abs_end,emi_start,emi_end;
                            r0=0.4,
                            theta=34.24,
                            #maxdt=20.,
                            scaletime=0.001,    
                            align=CAs)

# It is a good idea to save the data to a file:

file="tcf.dat"
println(" Writting data to file: ", file)
f = open(file,"w")
println(f," Time  2nd-Legendre TCF ")
for i in 1:length(t)
  println(f,t[i]," ",legendre[i]," ",tcf[i])
end
close(f)

println(" Plotting... ")
plot(t,legendre)
savefig("tcf.pdf")


