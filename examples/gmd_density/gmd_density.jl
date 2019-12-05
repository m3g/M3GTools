"""

Computing GMD contributions for the parts of the solute

"""

using M3GTools ; const m3g = M3GTools
using DelimitedFiles

#println(" Loading Plots... ")
#using Plots
#ENV["GKSwstype"]="nul" # This supresses the need of a display while plotting

mysim = m3g.init(psf="../gmd_solute/gmdfiles/structure.psf",
                 pdb="../gmd_solute/gmdfiles/structure.pdb",
                 vmd="vmd")

# Read GMD file corresponding to solute contributions
println(" Reading GMD solute contribution file... ")
gmd_solute_file = "../gmd_solute/gmdfiles/gmd-GMD_ATOM_SOLUTE_CONTRIB.dat"
gmd_solute = readdlm(gmd_solute_file,comments=true,comment_char='#')

# The distance is the first column of the gmd file, and the total gmd is the second column:
d = gmd_solute[:,1]
gmd_total = gmd_solute[:,2]

protein = m3g.select(mysim,"protein")

include("./Density.jl")
include("./setgrid.jl")

d = density3D(mysim,protein,gmd_solute_file)

## Of course, you can save the data for further analysis:
#
#output = open("./gmd_solute_contributions.dat","w")
#write(output,"# Distance Total Backbone Aliphatic Aromatic Polar\n")
#writedlm(output,zip(d,gmd_total,gmd_backbone,gmd_aliphatic,gmd_aromatic,gmd_polar))
#close(output)
#println(" Wrote file: ./gmd_solute_contribution.dat")
