
push!(LOAD_PATH,"../src")

using Namd
using Plots

function main()

  Namd.psf("./structure.psf")
  Namd.dcd("./structure.dcd")
  
  println("Reading masses...")
  mass = Namd.mass()

  println("Defining selections...")
  npopc, popc = Namd.select("resname POPC")
  nprot, prot = Namd.select("protein")
  
  println("Computing center of masses...")
  popc_cm = Array{Float32}(undef,Namd.nframes,3)
  prot_cm = Array{Float32}(undef,Namd.nframes,3)
  cm = Array{Float32}(undef,3)
  for i in 1:Namd.nframes
    sides, x, y, z = Namd.nextframe()
    cm = Namd.cm(popc,mass,x,y,z)
    for j in 1:3 
      popc_cm[i,j] = cm[j]
    end
    cm = Namd.cm(prot,mass,x,y,z)
    for j in 1:3 
      prot_cm[i,j] = cm[j]
    end
  end
  Namd.closedcd()

  println("Plotting...")
  diff = Vector{Float32}(undef,Namd.nframes)
  for i in 1:Namd.nframes
    diff[i] = prot_cm[i,3] - popc_cm[i,3]
  end
  x = [ i for i in 1:Namd.nframes ]
  plot(x,diff)
  savefig("example.pdf")

end; main()




