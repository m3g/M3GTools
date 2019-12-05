"""

3D represtation of the minimum distance distribution function

"""

using M3GTools ; const m3g = M3GTools

mysim = m3g.init(psf="../gmd_solute/gmdfiles/structure.psf",
                 pdb="../gmd_solute/gmdfiles/structure.pdb",
                 vmd="vmd")

# Selection that defines what the solute is

solute = m3g.select(mysim,"protein")

# File containting the contribution of the solute atoms to the gmd:

gmd_solute_file = "../gmd_solute/gmdfiles/gmd-GMD_ATOM_SOLUTE_CONTRIB.dat"

# Write PDB file containing the 3D respresentation of the GMD:

output = "./gmd3D.pdb"
gmd3D_write(mysim,solute,gmd_solute_file,output,)

