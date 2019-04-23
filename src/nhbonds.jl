#
# Count the number of hydrogen bonds between two selections
#

# If the coordinates of a PDB file read into the atoms vector is given
# Bonds are needed, because we need to know to which atoms the hydrogens
# are bound to
#

dotp(x,y) = x[1]*y[1]+x[2]*y[2]+x[3]*y[3]
vnorm(x) = sqrt(x[1]^2 + x[2]^2 + x[3]^2)
cosang(x,y) = dotp(x,y) / (vnorm(x)*vnorm(y))

function nhbonds(atoms,bonds,sel1,sel2; d = 3., angle = 20. )

  dhbond2 = d^2
  coshbond = cos(pi*angle/180.)

  n1 = length(sel1)
  n2 = length(sel2)

  nhbonds = 0
  for isel in 1:n1
    i = sel1[isel]
    for jsel in 1:n2
      j = sel2[jsel]
      # If one of the atoms is not electronegative, skip
      if ! ( atoms[i].element in [ "N", "O", "F" ] && atoms[j].element in [ "N", "O", "F" ] )
        continue
      end
      # If the distance is not compatible, skip
      d2 = (atoms[i].coor[1] - atoms[j].coor[1])^2 +
           (atoms[i].coor[2] - atoms[j].coor[2])^2 +
           (atoms[i].coor[3] - atoms[j].coor[3])^2
      if d2 > dhbond2
        continue
      end
      # Check if a hydrogen is bound to atom i
      foundhbond = false
      ibonds = filter( bond -> ( bond.i == i || bond.j == i ), bonds )
      for bond in ibonds
        if atoms[bond.i].element == "H"
          ib = bond.i
        elseif atoms[bond.j].element == "H" 
          ib = bond.j
        else  
          continue
        end
        v1H = [ atoms[ib].coor[1] - atoms[i].coor[1] ,
                atoms[ib].coor[2] - atoms[i].coor[2] ,
                atoms[ib].coor[3] - atoms[i].coor[3] ]
        vH2 = [ atoms[j].coor[1] - atoms[ib].coor[1] ,
                atoms[j].coor[2] - atoms[ib].coor[2] ,
                atoms[j].coor[3] - atoms[ib].coor[3] ]
        if cosang(v1H,vH2) > coshbond
          nhbonds = nhbonds + 1
          foundhbond = true
          break
        end
      end
      # If the hbond is already found between these to heavy atoms, skip next
      if foundhbond 
        continue
      end
      # Check if a hydrogen is bound to atom j
      jbonds = filter( bond -> ( bond.i == j || bond.j == j ), bonds )
      for bond in jbonds
        if atoms[bond.i].element == "H"
          ib = bond.i
        elseif atoms[bond.j].element == "H" 
          ib = bond.j
        else
          continue
        end
        v1H = [ atoms[ib].coor[1] - atoms[j].coor[1] ,
                atoms[ib].coor[2] - atoms[j].coor[2] ,
                atoms[ib].coor[3] - atoms[j].coor[3] ]
        vH2 = [ atoms[i].coor[1] - atoms[ib].coor[1] ,
                atoms[i].coor[2] - atoms[ib].coor[2] ,
                atoms[i].coor[3] - atoms[ib].coor[3] ]
        if cosang(v1H,vH2) > coshbond
          nhbonds = nhbonds + 1
          break
        end
      end
    end
  end

  return nhbonds
end

