"""

Computing GMD contributions for the parts of the solute

"""
function getgmd( gmd :: Array{}, sel1 :: String, sel2 :: String )

  # gmd is the array resulting from reading the gmd contributions per atom
  # for the solute, produced by the gmd.f90 program of mdanalysis

  # The first selection corresponds the the total solute, the second
  # selection to the section of the solute that will be considered here 

  sel1 = Namd.select(mysim,sel1)
  sel2 = Namd.select(mysim,sel2)

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

end function getgmd

getgmd( gmd :: Array{}; data = "none"; get = "none" ) = return getgmd( gmd, data, get )

