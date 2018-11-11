
push!(LOAD_PATH,"../src")

using Namd

println(" Loading Plots... ")
using Plots

function main()

  # Initialize simulation data and path to VMD executable
  println(" Reading simulation data ... ")
  Namd.init(psf="./structure.psf",
            dcd="./structure.dcd",
            vmd="/usr/bin/vmd")

  println(" Defining selections... ")
  popc = Namd.select("resname POPC")
  prot = Namd.select("protein")
  
  println(" Computing center of masses... ")
  popc_cm = Array{Float32}(undef,Namd.nframes,3)
  prot_cm = Array{Float32}(undef,Namd.nframes,3)
  cm = Array{Float32}(undef,3)
  for i in 1:Namd.nframes
    sides, x, y, z = Namd.nextframe()
    cm = Namd.cm(popc,Namd.mass,x,y,z)
    for j in 1:3 
      popc_cm[i,j] = cm[j]
    end
    cm = Namd.cm(prot,Namd.mass,x,y,z)
    for j in 1:3 
      prot_cm[i,j] = cm[j]
    end
  end
  Namd.closedcd()

  # Printing the difference in the Z-coordinate of the center of
  # masses of the selections

  println(" Plotting... ")
  diff = Vector{Float32}(undef,Namd.nframes)
  for i in 1:Namd.nframes
    diff[i] = prot_cm[i,3] - popc_cm[i,3]
  end
  x = [ i for i in 1:Namd.nframes ]
  plot(x,diff)
  savefig("example.pdf")

end; main()




