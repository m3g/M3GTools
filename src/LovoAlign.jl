
module LovoAlign

  export Pair

  #struct Pair
  #  p1 :: String
  #  p2 :: String
  #  TMScore :: Float64
  #end
  #Pair(;p1=p1,p2=p2) = Pair(p1,p2)

  function nonseq(pdb1,pdb2,sel1,sel2,output;
                  vmd_exec="vmd",
                  lovoalign_exec="lovoalign",
                  print=false)

    if pdb1 == output || pdb2 == output
      error("ERROR: Output file has the same name of one of input files.")
    end

    vmd_input = open("./LovoAlign_VMDINPUT_TMP.VMD","w")
    write(vmd_input,"""
                    mol new \"$pdb1\"
                    set all [ atomselect top all ]
                    \$all set beta 0
                    set sel [ atomselect top \"$sel1\" ]
                    \$sel set beta 1
                    \$all writepdb ./LovoAlign_sel1.pdb
                    \$all delete
                    \$sel delete
                    
                    mol new \"$pdb2\"
                    set all [ atomselect top all ]
                    \$all set beta 0
                    set sel [ atomselect top \"$sel2\" ]
                    \$sel set beta 1
                    \$all writepdb ./LovoAlign_sel2.pdb
                    \$all delete
                    \$sel delete
                    
                    exit
                    """)
    close(vmd_input)

    # Run VMD
    vmd_output = read(`$vmd_exec -dispdev text -e ./LovoAlign_VMDINPUT_TMP.VMD`,String)

    # Run LovoAlign
    lovoalign_output=read(
     `$lovoalign_exec -p1 LovoAlign_sel1.pdb -p2 LovoAlign_sel2.pdb -g 0. -beta1 -beta2 -m 3 -all -nglobal 50 -dtri 10. -o $output`,
     String)
 
    # Remove temporary files
    run(`\rm -f ./LovoAlign_VMDINPUT_TMP.VMD`)
    run(`\rm -f ./LovoAlign_sel1.pdb`)
    run(`\rm -f ./LovoAlign_sel2.pdb`)

    if print
      println(lovoalign_output)
      return
    else
      println("Wrote $output file with $pdb1 aligned to $pdb2")
      return
    end

  end

end

