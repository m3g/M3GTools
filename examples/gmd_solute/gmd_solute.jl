"""

Computing GMD contributions for the parts of the solute

"""

push!(LOAD_PATH,"../../")

using Namd
using DelimitedFiles

println(" Loading Plots... ")
using Plots
ENV["GKSwstype"]="nul" # This supresses the need of a display while plotting

# Data for my simulation (we only need the topoloy file here, might be PSF or GRO
topology = "./gmdfiles/structure.psf"

# Read GMD file corresponding to solute contributions
println(" Reading GMD solute contribution file... ")
gmd_solute = readdlm("./gmdfiles/gmd-GMD_ATOM_SOLUTE_CONTRIB.dat",comments=true,comment_char='#')

# The distance is the first column of the gmd file, and the total gmd is the second column:
d = gmd_solute[:,1]
gmd_total = gmd_solute[:,2]

# Define the solute selection used for computation of the GMD here:
gmd_backbone = Namd.gmdget(topology,gmd_solute,data="protein",get="protein and (not sidechain)")
gmd_aliphatic = Namd.gmdget(topology,gmd_solute,data="protein",get="protein and sidechain and aliphatic")
gmd_aromatic = Namd.gmdget(topology,gmd_solute,data="protein",get="protein and sidechain and aromatic")
gmd_polar = Namd.gmdget(topology,gmd_solute,data="protein",get="protein and sidechain and polar")

# Plot the results
plot(xlim=[0,8],xlabel="Distance / Angstrom",ylabel="MDDF")
plot!(d,gmd_total,label="Total",linewidth=2)
plot!(d,gmd_backbone,label="Backbone",linewidth=2)
plot!(d,gmd_aliphatic,label="Aliphatic",linewidth=2)
plot!(d,gmd_aromatic,label="Aromatic",linewidth=2)
plot!(d,gmd_polar,label="Polar",linewidth=2)
savefig("./gmd_solute.pdf")

# Of course, you can save the data for further analysis:

output = open("./gmd_solute_contributions.dat","w")
write(output,"# Distance Total Backbone Aliphatic Aromatic Polar\n")
writedlm(output,zip(d,gmd_total,gmd_backbone,gmd_aliphatic,gmd_aromatic,gmd_polar))
close(output)
println(" Wrote file: ./gmd_solute_contribution.dat")
