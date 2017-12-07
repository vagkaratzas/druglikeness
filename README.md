# druglikeness
This R script receives drug smiles per lines as an input and prints/returns the drug violations and calculated results, in two separate files, per drug/rule in a 2D matrix for each file.

Required libraries:
ChemmineR (https://www.bioconductor.org/packages/devel/bioc/vignettes/ChemmineR/inst/doc/ChemmineR.html#smiles-import)
ChemmineOB (https://www.bioconductor.org/packages/release/bioc/html/ChemmineOB.html)

Run example: Rscript druglikeness.R input-smiles.txt

Rscript must be on the environment path, else type whole path on command line while executing.
