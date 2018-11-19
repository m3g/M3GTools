
module LovoAlign

  using Random

  #export Pair
  #struct Pair
  #  p1 :: String
  #  p2 :: String
  #  TMScore :: Float64
  #end
  #Pair(;p1=p1,p2=p2) = Pair(p1,p2)

  #
  # Perform a alignment from coordinates, not pdb names
  #

  function align(x1,x2,y;
                 lovoalign_exec="lovoalign",debug=false)

   # Write temporary PDB file to lovoalign
 
   n1 = size(x1)[1]
   file = open("1.pdb","w")
   for i in 1:n1
     write(file,"ATOM 

   end
   
voltar

  end


  #
  # Perform a general, non-sequential alignment, using the Non-Bijective method
  #

  function nonseq(pdb1,pdb2,sel1,sel2,output;
                  vmd_exec="vmd",
                  lovoalign_exec="lovoalign",
                  debug=false)

    checkpdbs(pdb1,pdb2,output)

    # Define first selection using VMD
    sel1_pdb = set_selection_pdb(pdb1,sel1,vmd_exec=vmd_exec)

    # Define second selection using VMD
    sel2_pdb = set_selection_pdb(pdb2,sel2,vmd_exec=vmd_exec,debug=debug)

    # Run LovoAlign for general non-sequential alignment
    lovoalign_output=read(
     `$lovoalign_exec -p1 $sel1_pdb -p2 $sel2_pdb -g 0. -beta1 -beta2 -m 3 -all -nglobal 50 -dtri 10. -o $output`,
     String)
 
    # Remove temporary files
    run(`\rm -f $sel1_pdb`)
    run(`\rm -f $sel2_pdb`)

    if debug
      println(lovoalign_output)
      return
    else
      println("Wrote $output file with $pdb1 aligned to $pdb2")
      return
    end

  end

  #
  # Perform a maximization of the TM-score 
  #

  function TMScore(pdb1,pdb2,sel1,sel2,output;
                   vmd_exec="vmd",
                   lovoalign_exec="lovoalign",
                   debug=false)

    checkpdbs(pdb1,pdb2,output)

    # Define first selection using VMD
    sel1_pdb = set_selection_pdb(pdb1,sel1,vmd_exec=vmd_exec)

    # Define second selection using VMD
    sel2_pdb = set_selection_pdb(pdb2,sel2,vmd_exec=vmd_exec,debug=debug)

    # Run LovoAlign using default (TMscore) computation
    lovoalign_output=read(
     `$lovoalign_exec -p1 $sel1_pdb -p2 $sel2_pdb -o $output`,
     String)
 
    # Remove temporary files
    run(`\rm -f $sel1_pdb`)
    run(`\rm -f $sel2_pdb`)

    if debug
      println(lovoalign_output)
      return
    else
      println("Wrote $output file with $pdb1 aligned to $pdb2")
      return
    end

  end

  #
  # Function that checks some properteis of the PDB files
  #

  function checkpdbs(pdb1,pdb2,output)

    if pdb1 == output || pdb2 == output
      error("ERROR: Output file has the same name of one of input files.")
    end

  end

  #
  # Function that calls VMD to define the selected atoms
  #

  function set_selection_pdb(pdb,sel;vmd_exec="vmd",debug=false)

    rnd_chars=Random.randstring(8)
    filename = "./LovoAlign__$(rnd_chars)_$pdb"
    vmd_input_name = "./LovoAlign_$(rnd_chars).vmd"
    vmd_input = open(vmd_input_name,"w")
    write(vmd_input,"""
                    mol new \"$pdb\"
                    set all [ atomselect top all ]
                    \$all set beta 0
                    set sel [ atomselect top \"$sel\" ]
                    \$sel set beta 1
                    \$all writepdb $filename
                    \$all delete
                    \$sel delete
                    exit
                    """)
    close(vmd_input)
    vmd_output = read(`$vmd_exec -dispdev text -e $vmd_input_name`,String)
    if debug
      println(vmd_output)
    end
    run(`\rm -f $vmd_input_name`)
    return filename

  end

end

