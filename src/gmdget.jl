"""

Computing GMD contributions for the parts of the solute

"""

function gmdget(topology :: String, gmd :: Array{}, sel1 :: String, sel2 :: String; vmd="vmd" )

  # gmd is the array resulting from reading the gmd contributions per atom
  # for the solute, produced by the gmd.f90 program of mdanalysis

  # The first selection corresponds the the total solute, the second
  # selection to the section of the solute that will be considered here 

  sel1 = M3GTools.select(topology,sel1,vmd=vmd)
  sel2 = M3GTools.select(topology,sel2,vmd=vmd)

  # To do so, define a vector that contains the index in "solute" of the atoms
  # of the "backbone" selection. We sum 2 because this will correspond directly
  # to the columns of the gmd_solute array.

  iatoms = similar(sel2);
  for i in 1:length(sel2)
    iatoms[i] = findfirst(isequal(sel2[i]),sel1) + 2
  end
  
  # Sum all contributions of the gmd_solute array for iatoms:
  gmd_get = zeros(length(gmd[:,1]));
  for i in 1:length(gmd[:,1])
    for j in 1:length(sel2)
      gmd_get[i] = gmd_get[i] + gmd[i,iatoms[j]]
    end
  end

  return gmd_get

end

gmdget(topology :: String, gmd :: Array{}; data = "none", get = "none", vmd="vmd" ) = gmdget(topology, gmd, data, get, vmd=vmd)

gmdget(simulation :: Simulation, gmd :: Array{}, data, get) = gmdget(simulation, gmd, data, get, vmd=simulation.vmd)
gmdget(simulation :: Simulation, gmd :: Array{}; data = "none", get = "none" ) = gmdget(simulation.psf, gmd, data, get, vmd=simulation.vmd)

