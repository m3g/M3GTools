"""

Computing GMD contributions for the parts of the solute

"""

push!(LOAD_PATH,"/home/leandro/programs/namdjl")

using Namd
using DelimitedFiles

println(" Loading Plots... ")
using Plots
ENV["GKSwstype"]="nul" # This supresses the need of a display while plotting

# Data for my simulation (we only need the PSF file here)
mysim = Namd.init(psf="./teste.psf",vmd="vmd")

# Initialize simulation data and path to VMD executable
println(" Reading GMD solute contribution file... ")
gmd_solute = readdlm("./teste-GMD_ATOM_SOLUTE_CONTRIB.dat",comments=true,comment_char='#')

# The distance is the first column of the gmd file, and the total gmd is the second column:
d = gmd_solute[:,1]
gmd_total = gmd_solute[:,2]

# Define the solute selection used for computation of the GMD here:

gmd_backbone = Namd.gmdget(gmd_solute,data="protein",get="protein and backbone")
gmd_aliphatic = Namd.gmdget(gmd_solute,data="protein",get="protein and aliphatic")
gmd_charged = Namd.gmdget(gmd_solute,data="protein",get="protein and charged")

# Plot the results

plot(d,gmd_total,label="Total")
plot!(d,gmd_backbone,label="Backbone")
plot!(d,gmd_aliphatic,label="Aliphatic")
plot!(d,gmd_charged,label="Charged")
savefig("./gmd_solute.pdf")

# Of course, you can save the data for further analysis:

output = open("./gmd_solute_contributions.dat","w")
for i in 1:length(d)
  println(data,"$d[i] $gmd_total[i], $gmd_backbone[i], $gmd_aliphatic[i], $gmd_charged[i]")
end
close(output)


