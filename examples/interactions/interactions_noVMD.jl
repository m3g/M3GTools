#
# git clone https://github.com/mcubeg/namdjl
#

# In this example we compute the energy without using the complete Namd.jl package, and
# without calling VMD in background

# Path to PDBEnergy module
push!(LOAD_PATH, "../../");

# Load PDB energy module
using PDBEnergy

# Load PSF file data
atoms, bonds, angles, dihedrals, impropers = readpsf("../simulation_files/structure.psf");

# Load coordinates from PDB
getpdbcoords!("../simulation_files/structure.pdb",atoms);

# Select atoms of residue number 34
sel1 = [ atom.index for atom in filter( atom -> atom.residue == 34, atoms) ];

# Select atoms of residue number 30
sel2 = [ atom.index for atom in filter( atom -> atom.residue == 30, atoms) ];

# Compute electrostatic interaction between these two selections
elecenergy = coulomb(atoms,sel1,sel2);

# Read parameter files (which contain eps and sig vdW parameters)
parfiles = [ "/home/leandro/programs/toppar/charmm/par_all36_prot.prm", 
             "/home/leandro/programs/toppar/charmm/toppar_water_ions.str", 
             "/home/leandro/programs/toppar/charmm/par_all36_lipid.prm" ]
readprm!(parfiles,atoms)

# Compute vdw interaction between these two selections
vdwenergy = vdw(atoms,sel1,sel2);

println("Electrostatic interaction between residues 30 and 34 = $elecenergy")
println("vdW interaction between residues 30 and 34 = $vdwenergy") 

# Compute the number of hydrogen bonds between two selections
nhb = nhbonds(atoms,bonds,sel1,sel2)
println("Number of H-bonds between residues 30 and 34 = $nhb")


