# 
# Select atoms using vmd selection syntax, with vmd in background
#

function select(selection;update=false,
                          sides::Vector{Float64},
                          x::Vector{Float32},
                          y::Vector{Float32},
                          z::Vector{Float32})

  index_list = String
  readnext = false

  # If the selection has to be updated at every frame, we need
  # to write a temporary DCD file here to VMD

  if update
    writedcd(natoms,1,dcdaxis,sides,x,y,z;filename="Namdjl_DCDTEMP.dcd")
  end

  vmd_input = Base.open("./VMDINPUT_TMP.VMD","w")
  Base.write(vmd_input,"mol new \"$psffile\" \n")
  Base.write(vmd_input,"set sel [ atomselect top \"$selection\" ] \n")
  Base.write(vmd_input,"puts \"INDEXLIST\" \n")
  Base.write(vmd_input,"set indexes [ \$sel get index ] \n")
  Base.write(vmd_input,"puts \"ENDINDEXLIST\" \n")
  Base.write(vmd_input,"exit \n")
  Base.close(vmd_input)

  vmd_output = read(`$vmd -dispdev text -e ./VMDINPUT_TMP.VMD`, String)

  for line in split(vmd_output,"\n")
    if readnext
      if line == "ENDINDEXLIST" 
        error("ERROR: Selection '$selection' does not contain any atom")
      end 
      index_list = line
      break
    end
    if line == "INDEXLIST" 
      readnext = true
    end
  end
  index_split = split(index_list)
  nsel = length(index_split)
  selection_indexes = Vector{Int64}(undef,nsel) 
  for i in 1:nsel
    selection_indexes[i] = parse(Int64,index_split[i]) + 1
  end

  run(`\rm -f ./VMDINPUT_TMP.VMD`)

  if ! update 
    println(" Selection '$selection' contains ",nsel," atoms ")
  end
  return selection_indexes

end

