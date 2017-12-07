#!/usr/bin/env Rscript
#setwd("C:/Users/vagos/Desktop/Bio/Paper-CoDReS")
suppressMessages(library("ChemmineR"))
suppressMessages(library("ChemmineOB"))

args = commandArgs(trailingOnly = TRUE) #allow use of args
if (length(args) != 1) { #check if correct number of args
	stop("Need one input argument file containing drug SMILES.", call.=FALSE) #if not exactly one arg, throw error
}
smiles_matrix <- read.csv(args[1], header = FALSE) #smiles-input.txt
sdf <- apply(smiles_matrix, 1, smiles2sdf) #1 = rows, convert smiles to sdfs
violations <- (matrix(nrow = length(smiles_matrix[,1]))) #initialization of output violations' matrix
mol_weights <- lapply(sdf, MW) #apply on list, find molecular weights of sdfs
mol_weights_matrix <- matrix(unlist(mol_weights)) #stats on weights
mol_weights_violations <- lapply(mol_weights, function(x) {if(x>500) return(1) else return(0)})#(if(x > 500) return 1 else return 0))
mol_weights_mol_weights_violations_matrix <- matrix(unlist(mol_weights_violations))
violations <- cbind(smiles_matrix, mol_weights_mol_weights_violations_matrix)
stats <- cbind(smiles_matrix, mol_weights_matrix)

colnames(violations) <- c("SMILES", "Molecular Weights > 500")
colnames(stats) <- c("SMILES", "Molecular Weights")
print(violations)
print(stats)
write.table(violations, "druglikeness-violations.csv", row.names=FALSE, col.names=TRUE, sep=",")
write.table(stats, "druglikeness-stats.csv", row.names=FALSE, col.names=TRUE, sep=",")

