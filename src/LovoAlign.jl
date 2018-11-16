
module LovoAlign

  function nonseq(pdb1,pdb2,sel1,sel2,output;
                  vmd_exec="vmd",
                  lovoalign_exec="lovoalign")

    vmd_input = open("./LovoAlign_VMDINPUT_TMP.VMD","w")
    write("""
          mol new $pdb1
          set all [ atomselect top all ]
          $all set beta 0
          set sel [ atomselect top "$sel1" ]
          $sel set beta 1
          $all writepdb ./LovoAlign_sel1.pdb
          $all delete
          $sel delete
          
          mol new $pdb2
          set all [ atomselect top all ]
          $all set beta 0
          set sel [ atomselect top "$sel2" ]
          $sel set beta 1
          $all writepdb ./LovoAlign_sel2.pdb
          $all delete
          $sel delete
          
          exit
          """)
    close(vmd_input)

    # Run VMD
    run(`$vmd_exec -dispdev text -e ./LovoAlign_VMDINPUT_TMP.VMD`)

    # Run LovoAlign
    run(`$lovoalign_exec -p1 LovoAlign_sel1.pdb -p2 LovoAlign_sel2.pdb -g 0. -beta1 -beta2 -m 3 -all -nglobal 10 -o $output`)

    # Remove temporary files
    run(`\rm -f ./LovoAlign_VMDINPUT_TMP.VMD`)
    run(`\rm -f ./LovoAlign_sel1.pdb`)
    run(`\rm -f ./LovoAlign_sel2.pdb`)

  end

end

