## GMD 3D: Writes a three-dimensional (3D) representation of the minimum-distance distribution

In this example we show how to represent the MDDF in 3D, to visualize the
distribution on top of the solute structure.

The MDDF must be computed with the GMD module of the MDAnalysis suite of
packages (version 19.148 or greater), which is available here:

http://leandro.iqm.unicamp.br/mdanalysis/gmd

In this example, we had protein in a solution containing the stabilizer
TMAO, and we computed the MDDF between the protein and TMAO (for more
details, see: http://pubs.acs.org/doi/abs/10.1021/acs.jctc.7b00599).

To generate this data, just follow the example provided in
<a href=https://github.com/mcubeg/namdjl/blob/master/examples/gmd3D/gmd3D.jl>`gmd3D.jl`</a>. 
You will need the file created by the GMD module of
MDAnalysis mentioned above, containing the contribution of each solute
atom the distribution function. This file, in this example, is called:
``` 
gmd-GMD_ATOM_SOLUTE_CONTRIB.dat
```

First, we need to read the simulation PDB and PSF files, using:
```
mysim = m3g.init(psf="../gmd_solute/gmdfiles/structure.psf",
                 pdb="../gmd_solute/gmdfiles/structure.pdb",
                 vmd="vmd")
```                   
where `vmd` the path to the VMD executable.

In this case, the solute is the protein in the system:

```                   
solute = m3g.select(mysim,"protein")
```                   

This is the name of the file containing the contributions of the solute atoms
to the MDDF:
```                   
gmd_solute_file = "../gmd_solute/gmdfiles/gmd-GMD_ATOM_SOLUTE_CONTRIB.dat"
```                   

Given the above data, just call this function to obtain a 3D representation
of the MDDF, which might be loaded with the protein structure
in any visualization package: 

```                   
output = "./gmd3D.pdb"
gmd3D_write(mysim,solute,gmd_solute_file,output,)
```                   

The PDB file generated has the following information, for example:

```                   
ATOM      1  HZ3 LYS    25     -20.680   1.439  10.131  6.93 10.97
ATOM      2  HZ3 LYS    25     -20.680   1.439  10.631  6.93  9.57
ATOM      3  HZ3 LYS    25     -20.680   1.439  11.131  6.95  5.15
ATOM      4  HZ3 LYS    25     -20.680   1.939  12.131  6.96  2.11
```                   

The positions correspond to the positions in the 3D grid, around the solute, where 
the MDDF is being reported. The atom, in this case `HZ3 LYS 25` is the closest atom
of the protein to that point in space. That is, that is the atom which will count
in the MDDF for that position in space. Thus, if the minimum distance between
the solute and the solvent molecule occurs at that point in space, that would
be the atom of the solute considered.     

The occupancy field (with `6.93`, for example) contains the distance of that grid point from
the closest solute atom. 

Finally, the B-factor field (with `10.97`, for example) contains the contribution
of that protein atom to the total MDDF at that distance, normalized such that the 
maximum contribution is 100.   

Essentially, this PDB file will contain atoms at the positions around the protein that contribute
to the MDDF at each distance. For example, if we plot the points for which the occupancy is 
less than 2.0, we get the positions around the protein where the solvent molecules
formed hydrogen bonds: 

<p align="center" width=100%>
<img width=30% src="https://github.com/mcubeg/M3GTools/raw/master/examples/gmd3D/figures/occup_lt_2.png">
<img width=30% src="https://github.com/mcubeg/M3GTools/raw/master/examples/gmd3D/figures/occup_lt_2_atomtype.png">
<img width=30% src="https://github.com/mcubeg/M3GTools/raw/master/examples/gmd3D/figures/occup_lt_2_resname.png">
</p>

(Figures made with VMD, using the "QuickSurf" representation for the grid). 

In the figure of the left above, from blue to red see the position that contribute more to the hydrogen-bonding
peak of the MDDF, because each point is colored by its B-factor. In the
figure of the center we show which type of protein atom is involved in
the interactions (red: Oxygen, white: Hydrogen), and thus the
acceptor/donnor character of the hydrogen bonds becomes clear at each
position. Finally, in the figure of the right we paint the grid by the
residue type involved in the interactions (for example: Yellow: Serine
and Red: Aspartic acid). 










