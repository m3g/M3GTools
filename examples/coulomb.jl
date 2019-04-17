
# Path to PDBEnergy module
push!(LOAD_PATH, "../");

# Load PDB energy module
using PDBEnergy

# Load PSF file data
atoms, bonds, angles, dihedrals, impropers = readpsf("./structure.psf")

# Load coordinates from PDB
getpdbcoords!("./structure.pdb",atoms)

# Select atoms of residue number 17
sel1 = [ atom.index for atom in filter( atom -> atom.residue == 17, atoms) ];

# Select atoms of residue number 32
sel2 = [ atom.index for atom in filter( atom -> atom.residue == 19, atoms) ];

# Compute electrostatic interaction between these two selections
q = coulomb(atoms,sel1,sel2)

println("Electrostatic interaction between residues 17 and 19 = $q")

