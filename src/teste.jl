using FortranFiles
include("./structures.jl")
include("./readpsf.jl")
include("./readprm.jl")

atoms, bonds, angles, dihedrals, impropers = readpsf("../examples/structure.psf")

parfiles = [ "/home/leandro/programs/toppar/charmm/par_all36_prot.prm",
             "/home/leandro/programs/toppar/charmm/toppar_water_ions.str",
             "/home/leandro/programs/toppar/charmm/par_all36_lipid.prm" ]

readprm!(parfiles,atoms,angles)
