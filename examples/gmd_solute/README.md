## GMD_SOLUTE: Extract the contributions of parts of the solute to the minimum-distance distribution function

This example shows how to extract from the total minimum-distance
distribution function (MDDF) the contribution of specific parts of the
solute. 

For example, we will decompose the MDDF of a protein-solvent
distribution function into the contributions of the backbone and some
specific side-chains. 

The MDDF must be computed with the GMD module of the MDAnalysis suite of
packages (version 19.148 or greater), which is available here:

http://leandro.iqm.unicamp.br/mdanalysis/gmd

In this example, we had protein in a solution containing the stabilzer
TMAO, and we computed the MDDF between the protein and TMAO (for more
details, see: http://pubs.acs.org/doi/abs/10.1021/acs.jctc.7b00599). 

The expected result from this analysis is the decomposition of the MDDF
into different subsets of the solute structure. Here, we decompose the
Protein-TMAO MDDF into the contributions of the protein backbone,
aliphatic, aromatic and polar side chains. 

<p align="center">
<img src="https://github.com/mcubeg/namdjl/raw/master/examples/gmd_solute/gmdfiles/gmd_total.png">
</p>

To generate this data, just follow the example provided in
<a href=https://github.com/mcubeg/namdjl/blob/master/examples/gmd_solute/gmd_solute.jl>`gmd_solute.jl`</a>. You will need the file created by the GMD module of
MDAnalysis mentioned above, containing the contribution of each solute
atom the distribution function. This file, in this example, is called:
``` 
gmd-GMD_ATOM_SOLUTE_CONTRIB.dat
```
This file is loaded with the command (which requires the `DelimitedFiles` package):

```
gmd_solute = readdlm("./gmdfiles/gmd-GMD_ATOM_SOLUTE_CONTRIB.dat",comments=true,comment_char='#')
```
The first column of the data read is the distance, the second column is
the total MDDF, and the other columns are the contributions of each atom
of the solute to the MDDF. We can get the first two columns with: 
```
d = gmd_solute[:,1]
gmd_total = gmd_solute[:,2]
```
To extract the contribution of a specific set of atoms given by a
selection with VMD syntax, use the `gmdget` function, as in the example:
```
topology = "/structure.psf"
gmd_backbone = Namd.gmdget(topology,gmd_solute,data="protein",get="protein and (not sidechain)")
```
where `topology` is the topology file for the system, which might be of
PSF, GRO or any other format which VMD understands.

Given the desired selections, a plot similar to the above can be
produced, with, for example,
```
plot(xlim=[0,8],xlabel="Distance / Angstrom",ylabel="MDDF")
plot!(d,gmd_total,label="Total",linewidth=2)
plot!(d,gmd_backbone,label="Backbone",linewidth=2)
savefig("./gmd_solute.pdf")
```
In which the contribution of the interactions of the solvent with each
part of the solute can be clearly discriminated.





