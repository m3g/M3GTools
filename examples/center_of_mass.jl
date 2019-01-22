"""

In this example, the centers of mass of the protein and of the POPC residues
are computed, and the difference between the z-coordinates of them is ploted.

"""

push!(LOAD_PATH,"../src")

using Namd

println(" Loading Plots... ")
using Plots

function main()

  # Initialize simulation data and path to VMD executable

  println(" Reading simulation data ... ")

  mysim = Namd.init(psf="./structure.psf",
                    dcd="./structure.dcd",
                    vmd="vmd")

  println(" Defining selections... ")
  popc = Namd.select(mysim,"resname POPC")
  prot = Namd.select(mysim,"protein")

  println(" Computing center of masses... ")
  popc_cm = Array{Float32}(undef,mysim.nframes,3)
  prot_cm = Array{Float32}(undef,mysim.nframes,3)
  cm = Array{Float32}(undef,3)

  for i in 1:mysim.nframes
    sides, x, y, z = Namd.nextframe(mysim)

    cm = Namd.cm(popc,mysim.mass,x,y,z)
    for j in 1:3 
      popc_cm[i,j] = cm[j]
    end
    cm = Namd.cm(prot,mysim.mass,x,y,z)
    for j in 1:3 
      prot_cm[i,j] = cm[j]
    end

  end
  Namd.closedcd(mysim)

  # Printing the difference in the Z-coordinate of the center of
  # masses of the selections

  println(" Plotting... ")
  diff = Vector{Float32}(undef,mysim.nframes)
  for i in 1:mysim.nframes
    diff[i] = prot_cm[i,3] - popc_cm[i,3]
  end
  x = [ i for i in 1:mysim.nframes ]
  plot(x,diff)
  savefig("example1.pdf")

end; main()




