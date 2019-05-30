## GMD_SOLUTE: Extract the contributions of parts of the solute to the minimum-distance distribution function

This example shows how to extract from the total minimum-distance distribution function (MDDF) the contribution of specific 
parts of the solute. 

For example, we will decompose the MDDF of a protein-solvent distribution function into the contributions of the backbone
and some specific side-chains. 

The MDDF must be computed with the GMD module of the MDAnalysis suite of packages, which is available here:

http://leandro.iqm.unicamp.br/mdanalysis/gmd

In this example, we had protein in a solution containing the stabilzer TMAO, and we computed the MDDF between the protein
and TMAO (for more details, see: http://pubs.acs.org/doi/abs/10.1021/acs.jctc.7b00599). 

The Protein-TMAO MDDF is shown below, and indicates that TMAO interacts with the protein surface:

<p align="center">
<img src="https://github.com/mcubeg/namdjl/blob/master/examples/gmd_solute/gmdfiles/gmdtotal.png?raw=true">
</p>

