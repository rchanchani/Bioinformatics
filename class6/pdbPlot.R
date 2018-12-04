## Description: Reads a PDB file and from the chain and atom of interest,
##              extracts the beta factor and plots against residue if plot is
##              required.
## Usage: Provide PDB file of interest, the chain and atoms of interest, if plot is created (optional),
##        color of plot (optional), and whether or not bar graph of chain is included with line plot (optional)
##        Beta factor extracted from specified secondary structure and plotted (optional) against residues.
## Inputs: x - character - a PDB file of interest
##         chain - character - the chain of interest within the structure - default chain A
##         elety - character - a character vector of atom names you want to extract from chain;
##                            default "CA" central carbon.
##         plot - logical - whether or not a plot is generated
##         bar - logical - whether or not the bars representing chains are shown in the plot
##         col - character - color of the line plot
## Ouput: a line plot of the beta factor versus residues, if user selects TRUE for plot argument
pdb.plot <- function(x, chain = "A", elety = "CA", plot = FALSE, bar = TRUE, col = "black") {
  # loads the bio3d library to use read.pdb, plotb3
  library(bio3d)
  
  # read in pdb file; gives warning if file is already downloaded
  s <- read.pdb(x)
  
  # extract subset of atoms in a PDB file on the desired chain
  s.chain <- trim.pdb(s, chain = chain, elety = elety)
  
  # extract bata factor from the chain of interest
  s.b <- s.chain$atom$b
  
  # plot the beta factor versus the residue in a line plot with specified (optional) color
  if(plot) {
    
    # user specifies inclusion of secondary structure (chain) annotation
    if (bar) {
      return(plotb3(s.b, sse = s.chain, typ = "l", ylab = "Beta Factor", col = col,
                    main = "Residues vs. Beta Factor"))
    }
    
    # (default) removal of bar plot of secondary structure (chain) annotation
    else {
      return(plotb3(s.b, typ = "l", ylab = "Beta Factor", col = col,
                    main = "Residues vs. Beta Factor"))
    }
  }
}