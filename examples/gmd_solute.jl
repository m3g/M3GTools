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

# Define the solute selection used for computation of the GMD here:
solute = Namd.select(mysim,"protein")

# We will plot the contribution of the backbone for the GMD:
backbone = Namd.select(mysim,"protein and backbone")

# We will plot the contribution of the backbone for the GMD:
backbone = Namd.select(mysim,"protein and backbone")

# To do so, define a vector that contains the index in "solute" of the atoms
# of the "backbone" selection. We sum 2 because this will correspond directly
# to the columns of the gmd_solute array.
iatoms = similar(backbone);
for i in 1:length(backbone)
  iatoms[i] = findfirst(isequal(backbone[i]),solute) + 2
end

# Sum all contributions of the gmd_solute array for iatoms:
gmd_backbone = zeros(length(gmd_solute[:,1]));
for i in 1:length(gmd_solute[:,1])
  for j in 1:length(backbone)
    gmd_backbone[i] = gmd_backbone[i] + gmd_solute[i,iatoms[j]]
  end
end

# We will plot the contribution of the side chains for the GMD:
aliphatic = Namd.select(mysim,"protein and aliphatic")
iatoms = similar(sidechains);
for i in 1:length(aliphatic)
  iatoms[i] = findfirst(isequal(aliphatic[i]),solute) + 2
end
gmd_aliphatic = zeros(length(gmd_solute[:,1]));
for i in 1:length(gmd_solute[:,1])
  for j in 1:length(aliphatic)
    gmd_aliphatic[i] = gmd_aliphatic[i] + gmd_solute[i,iatoms[j]]
  end
end

# We will plot the contribution of the charged for the GMD:
charged = Namd.select(mysim,"protein and charged")

# To do so, define a vector that contains the index in "solute" of the atoms
# of the "backbone" selection. We sum 2 because this will correspond directly
# to the columns of the gmd_solute array.
iatoms = similar(charged);
for i in 1:length(charged)
  iatoms[i] = findfirst(isequal(charged[i]),solute) + 2
end

# Sum all contributions of the gmd_solute array for iatoms:
gmd_charged = zeros(length(gmd_solute[:,1]));
for i in 1:length(gmd_solute[:,1])
  for j in 1:length(charged)
    gmd_charged[i] = gmd_charged[i] + gmd_solute[i,iatoms[j]]
  end
end
plot(gmd_solute[:,1],gmd_solute[:,2],label="Total")
plot!(gmd_solute[:,1],gmd_backbone,label="Backbone")
plot!(gmd_solute[:,1],gmd_aliphatic,label="Aliphatic")
plot!(gmd_solute[:,1],gmd_charged,label="Charged")

savefig("./teste.pdf")






