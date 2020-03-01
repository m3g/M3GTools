#
# git clone https://github.com/m3g/M3GTools
#

# Load M3GTools
using M3GTools

# Load simulation data
mysim = M3GTools.init(psf="../simulation_files/structure.psf",
                  pdb="../simulation_files/structure.pdb",
                  parfiles = [ "/home/leandro/programs/toppar/charmm/par_all36_prot.prm", 
                               "/home/leandro/programs/toppar/charmm/toppar_water_ions.str", 
                               "/home/leandro/programs/toppar/charmm/par_all36_lipid.prm" ],
                  vmd="vmd")

# Select atoms of residue number 33 (be careful with the fact that VMD starts counting from zero):
sel1 = M3GTools.select(mysim,"protein and residue 33")

# Select atoms of residue number 29
sel2 = M3GTools.select(mysim,"protein and residue 29")

# Compute electrostatic interaction between these two selections
elecenergy = M3GTools.coulomb(mysim.atoms,sel1,sel2);

# Compute vdw interaction between these two selections
vdwenergy = M3GTools.vdw(mysim.atoms,sel1,sel2);

println("Electrostatic interaction between residues 30 and 34 = $elecenergy")
println("vdW interaction between residues 30 and 34 = $vdwenergy") 

# Compute the number of hydrogen bonds between two selections
nhb = M3GTools.nhbonds(mysim.atoms,mysim.bonds,sel1,sel2)
println("Number of H-bonds between residues 30 and 34 = $nhb")


