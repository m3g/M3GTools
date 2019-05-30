
"""

Computing GMD contributions for the parts of the solute

"""

push!(LOAD_PATH,"../../../")

using Namd
using DelimitedFiles
using LaTeXStrings

println(" Loading Plots... ")
using Plots
ENV["GKSwstype"]="nul" # This supresses the need of a display while plotting

# Data for my simulation (we only need the topoloy file here, might be PSF or GRO
topology = "./structure.psf"

# Read GMD file corresponding to solute contributions
println(" Reading GMD solute contribution file... ")
gmd_solute = readdlm("./gmd-GMD_ATOM_SOLUTE_CONTRIB.dat",comments=true,comment_char='#')

# The distance is the first column of the gmd file, and the total gmd is the second column:
d = gmd_solute[:,1]
gmd_total = gmd_solute[:,2]

# Define the solute selection used for computation of the GMD here:
gmd_backbone = Namd.gmdget(topology,gmd_solute,data="protein",get="protein and (backbone or name HN)")
gmd_aliphatic = Namd.gmdget(topology,gmd_solute,data="protein",get="protein and aliphatic")
gmd_charged = Namd.gmdget(topology,gmd_solute,data="protein",get="protein and charged")

# Plot the results
plot(d,gmd_total,label="Total",xlim=[0,8],xlabel=L"\textrm{Distance / \AA}",ylabel=L"\textrm{MDDF}")
plot!(d,gmd_backbone,label="Backbone")
plot!(d,gmd_aliphatic,label="Aliphatic")
plot!(d,gmd_charged,label="Charged")
savefig("./gmd_total.png")
println(" Created plot: ./gmd_total.png")

