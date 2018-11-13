Structural Bioinformatics II
================
Raghav Chanchani
11/13/2018

HIV Structural Bioinformatics
=============================

Splitting up PDB files into proteins and ligands
------------------------------------------------

We will first download and read in the 1hsg protein from PDB.

``` r
library(bio3d)
prot <- get.pdb("1hsg")
```

Then, we will read in the strucutre of the protein and the ligand bound to it, and split up these two components into separate files to attempt docking of ligands into the protein. This is done after removing the water molecules in the structure.

``` r
hiv <- read.pdb(prot)
prot <- trim.pdb(hiv, "protein")
lig <-trim.pdb(hiv, "protein")
write.pdb(prot, file="1hsg_protein.pdb")
write.pdb(lig, file="1hsg_ligand.pdb")
```