# 
# Select atoms using vmd selection syntax, with vmd in background
#

function select(simulation :: Simulation, selection :: String)

  local index_list :: String 
  local readnext :: Bool = false

  vmd_input = Base.open("./VMDINPUT_TMP.VMD","w")
  Base.write(vmd_input,"mol new \"$(simulation.psf)\" \n")
  Base.write(vmd_input,"set sel [ atomselect top \"$selection\" ] \n")
  Base.write(vmd_input,"puts \"INDEXLIST\" \n")
  Base.write(vmd_input,"set indexes [ \$sel get index ] \n")
  Base.write(vmd_input,"puts \"ENDINDEXLIST\" \n")
  Base.write(vmd_input,"exit \n")
  Base.close(vmd_input)

  vmd_output = read(`$(simulation.vmd) -dispdev text -e ./VMDINPUT_TMP.VMD`, String)

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

  return selection_indexes

end

