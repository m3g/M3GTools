
using Plots
using DelimitedFiles

data = readdlm("./gr.dat")

r = data[:,1]
g = data[:,2]

plot(r,g,label="",linewidth=2)
plot!(xlabel="r / Angstroms",ylabel="MDDF")

plot!(size=(300,300))

savefig("./gr.pdf")

