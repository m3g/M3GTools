"""

In this example, the distance between the C-terminal carbon
and a sodium atom is computed and ploted as a function of time

"""

using M3GTools

println(" Loading Plots... ")
using Plots
ENV["GKSwstype"]="nul" # This supresses the need of a display while plotting


# Initialize simulation data and path to VMD executable
println(" Reading simulation data ... ")
mysim = M3GTools.init(psf="../simulation_files/structure.psf",
                      dcd="../simulation_files/structure.dcd",
                      vmd="vmd")

println(" Defining selections... ")
cterm = M3GTools.select(mysim,"protein and resid 253 and name C")
sod = M3GTools.select(mysim,"index 20455")

# Vector that will contain the distances for each frame
d = Vector{Float64}(undef,mysim.nframes)

println(" Computing distances ... ")
for i in 1:mysim.nframes
  sides, x, y, z = M3GTools.nextframe(mysim)

  # Wrap coordinates of the sodium around cterm

  x_cterm = [ x[cterm[1]], y[cterm[1]], z[cterm[1]] ]
  M3GTools.wrap!(sides,x,y,z;center=x_cterm,sel=sod)

  # Compute distance

  d[i] = sqrt( ( x[cterm[1]] - x[sod[1]] )^2 + 
               ( y[cterm[1]] - y[sod[1]] )^2 + 
               ( z[cterm[1]] - z[sod[1]] )^2 ) 

end
M3GTools.closedcd(mysim)

println(" Plotting... ")

x = [ i for i in 1:mysim.nframes ]
plot(x,d,xlabel="Frame",ylabel="Distance")
savefig("distance.pdf")
println(" Created figure: distance.pdf")



