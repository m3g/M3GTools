"""

Computing GMD contributions for the parts of the solute

"""

push!(LOAD_PATH,"../")

using Namd
using DelimitedFiles

println(" Loading Plots... ")
using Plots
ENV["GKSwstype"]="nul" # This supresses the need of a display while plotting

# Data for my simulation (we only need the PSF file here)
mysim = Namd.init(psf="./gmdfiles/test.psf",vmd="vmd")

# Read GMD file corresponding to solute contributions
println(" Reading GMD solute contribution file... ")
gmd_solute = readdlm("./gmdfiles/gmd-GMD_ATOM_SOLUTE_CONTRIB.dat",comments=true,comment_char='#')

# The distance is the first column of the gmd file, and the total gmd is the second column:
d = gmd_solute[:,1]
gmd_total = gmd_solute[:,2]

# Define the solute selection used for computation of the GMD here:
gmd_backbone = Namd.gmdget(mysim,gmd_solute,data="protein",get="protein and backbone")
gmd_aliphatic = Namd.gmdget(mysim,gmd_solute,data="protein",get="protein and aliphatic")
gmd_charged = Namd.gmdget(mysim,gmd_solute,data="protein",get="protein and charged")

# Plot the results
plot(d,gmd_total,label="Total")
plot!(d,gmd_backbone,label="Backbone")
plot!(d,gmd_aliphatic,label="Aliphatic")
plot!(d,gmd_charged,label="Charged")
savefig("./gmd_solute.pdf")

# Of course, you can save the data for further analysis:

output = open("./gmd_solute_contributions.dat","w")
write(output,"# Distance Total Backbone Aliphatic Charged \n")
writedlm(output,zip(d,gmd_total,gmd_backbone,gmd_aliphatic,gmd_charged))
close(output)


