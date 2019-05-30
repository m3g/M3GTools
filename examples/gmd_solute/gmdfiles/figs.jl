
"""

Computing GMD contributions for the parts of the solute

"""

push!(LOAD_PATH,"../../../")

using Namd
using DelimitedFiles
using LaTeXStrings

function smooth(data)
  # smooth data
  n = length(data)
  y = copy(data)
  nav = 10
  for i in nav+1:n-nav
    av = 0.
    for j in i-nav:i+nav
      av = av + data[j]
    end
    y[i] = av / (2*nav)
  end
  for i in 1:nav
    y[i] = data[i]
    y[n-i] = data[n-i]
  end
  return y
end 

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

gmd_total = smooth(gmd_total)

# Define the solute selection used for computation of the GMD here:
gmd_backbone = Namd.gmdget(topology,gmd_solute,data="protein",get="protein and (not sidechain)")
gmd_aliphatic = Namd.gmdget(topology,gmd_solute,data="protein",get="protein and sidechain and aliphatic")
gmd_aromatic = Namd.gmdget(topology,gmd_solute,data="protein",get="protein and sidechain and aromatic")
gmd_polar = Namd.gmdget(topology,gmd_solute,data="protein",get="protein and sidechain and polar")

gmd_backbone = smooth(gmd_backbone)
gmd_aliphatic = smooth(gmd_aliphatic)
gmd_aromatic = smooth(gmd_aromatic)
gmd_polar = smooth(gmd_polar)

# Plot the results
plot(d,gmd_total,label="Total",xlim=[0,8],xlabel=L"\textrm{Distance / \AA}",ylabel=L"\textrm{MDDF}",linewidth=2)
plot!(d,gmd_backbone,label="Backbone",linewidth=2)
plot!(d,gmd_aliphatic,label="Aliphatic",linewidth=2)
plot!(d,gmd_aromatic,label="Aromatic",linewidth=2)
plot!(d,gmd_polar,label="Polar",linewidth=2)
plot!(fontsize=16)
savefig("./gmd_total.png")
println(" Created plot: ./gmd_total.png")

