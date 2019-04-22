#
# Count the number of hydrogen bonds between two selections
#

# If the coordinates of a PDB file read into the atoms vector is given
# Bonds are needed, because we need to know to which atoms the hydrogens
# are bound to

function hbonds(atoms,bonds,sel1,sel2; d = 3., angle = 20. )

  dhbond2 = d^2
  cosangle = cos(pi*angle/180.)

  n1 = length(sel1)
  n2 = length(sel2)
  boundtoH1 = zeros(Bool,n1)
  boundtoH2 = zeros(Bool,n2)

  nof1 = zeros(Bool,n1)
  for i in 1:n1
    if element(atoms[i].name) in [ "N", "O", "F" ]
      nof1[i] = true
      bonds = filter( bond -> ( bond.i == sel1[i] || bond.j == sel1[i] ), bonds )
      for bond in bonds
        if element(atoms[bond.j]) == "H") || elment(atoms[bond.i] == "H")
          boundtoH1[i] = true
        end
      end
    end
  end
  nof2 = zeros(Bool,n2)
  for i in 1:n2
    if element(atoms[i].name) in [ "N", "O", "F" ]
      nof2[i] = true
      bonds = filter( bond -> ( bond.i == sel2[i] || bond.j == sel2[i] ), bonds )
      for bond in bonds
        if element(atoms[bond.j]) == "H") || elment(atoms[bond.i] == "H")
          boundtoH1[i] = true
        end
      end
    end
  end

  for isel1 in 1:n1
    i = sel1[isel1]
    for isel2 in 1:n2
      j = sel2[isel2]
      if nof1[isel1] && nof2[isel2]
        if boundtoH1[isel1] || boundtoH2[isel2]

          # Compute the geometrical parameters

          d = atom[i].coor[1]-atom[j].coor[1])^2 +
              atom[i].coor[2]-atom[j].coor[2])^2 +
              atom[i].coor[3]-atom[j].coor[3])^2 
          if d < dhbond2
            # falta saber se o H está no meio! voltar

          end
        end 
      end
    end
  end

end

# If the coordinates of a simulation are given (this must be done using PBC
# and linked cells, to be efficient)
# ...


function element(type :: String)
  try
    i = parse(Int,type[1:1])
    return type[2:2]
  catch
    return type[1:1]
  end
end
