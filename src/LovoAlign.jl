
module LovoAlign

  using Random
  using PDB

  #export Pair
  #struct Pair
  #  p1 :: String
  #  p2 :: String
  #  TMScore :: Float64
  #end
  #Pair(;p1=p1,p2=p2) = Pair(p1,p2)

  #
  # Performs a alignment from coordinates, not pdb names
  #

  function align(x1,x2;
                 method=2,
                 lovoalign_exec="lovoalign",debug=false)

    rnd_chars=Random.randstring(8)

    # Write temporary PDB file to lovoalign
 
    n1 = size(x1)[1]
    file1_name = "LovoAlign_1_$(rnd_chars).pdb"
    file = open(file1_name,"w")
    for i in 1:n1
      line = PDB.printatom(i,"CA","GLY","A",i,x1[i,1],x1[i,2],x1[i,3])
      write(file,line,"\n")
    end
    close(file)
    
    n2 = size(x2)[1]
    file2_name = "LovoAlign_2_$(rnd_chars).pdb"
    file = open(file2_name,"w")
    for i in 1:n2
      line = PDB.printatom(i,"CA","GLY","A",i,x2[i,1],x2[i,2],x2[i,3])
      write(file,line,"\n")
    end
    close(file)

    if n1 < 15 || n2 < 15 
      method = 3
    end

    output="LovoAlign_3_$(rnd_chars).pdb"
    lovoalign_output=read(
     `$lovoalign_exec -p1 $file1_name -p2 $file2_name -m $method -o $output`,
     String)

    if debug
      println(lovoalign_output)
    end
 
    y = Array{Float32}(undef,n1,3)
    file = open(output,"r")
    i = 0
    for line in eachline(file)
      i = i + 1
      y[i,1] = parse(Float64,line[31:38])
      y[i,2] = parse(Float64,line[39:46])
      y[i,3] = parse(Float64,line[47:54])
    end
    close(output)

   run(`\rm -f $file1_name $file2_name $output`)

   return y
   
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

