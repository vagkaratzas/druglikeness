#!/usr/bin/env Rscript
#setwd("C:/Users/vagos/Desktop/Bio/Paper-CoDReS")
suppressMessages(library("ChemmineR"))
suppressMessages(library("ChemmineOB"))

args = commandArgs(trailingOnly = TRUE) #allow use of args
if (length(args) != 1) { #check if correct number of args
	stop("Need one input argument file containing drug SMILES.", call.=FALSE) #if not exactly one arg, throw error
}
smiles_matrix <- read.csv(args[1], header = FALSE) #smiles-input.txt
print("### Lipinski Rule of Five - Druglikeness Violations ###")
sdf <- apply(smiles_matrix, 1, smiles2sdf) #1 = rows, convert smiles to sdfs
openbabel <- lapply(sdf, propOB) #returns all the data we need to figure Lipinski's rules
hba1 <- matrix(0L, nrow = length(smiles_matrix[,1]))
hbd <- hba1 #initialization
logp <- hba1 #initialization
mw <- hba1 #initialization
for (i in 1:length(openbabel)){
	ob <- unlist(openbabel[i]) #in order to parse results
	hba1[i] <- ob[6] #hbond acceptors
	hbd[i] <- ob[8] #hbond donors
	logp[i] <- ob[9] #partition coefficient
	mw[i] <- ob[11] #molecular weight
}
stats <- cbind(smiles_matrix, hba1, hbd, logp, mw)
colnames(stats) <- c("SMILES", "Hydrogen Bond Donors", "Hydrogen Bond Acceptors", "Octanol/Water Partition Coefficient (logP)", "Molecular Weight")
hba1_violation <- lapply(hba1, function(x) {if(x>5) return(1) else return(0)})#figure if rule violated
hba1_violation_matrix <- matrix(unlist(hba1_violation)) #making amtrix from list to cbind later
hbd_violation <- lapply(hbd, function(x) {if(x>10) return(1) else return(0)})#figure if rule violated
hbd_violation_matrix <- matrix(unlist(hbd_violation)) #making amtrix from list to cbind later
logp_violation <- lapply(logp, function(x) {if(x>5) return(1) else return(0)})#figure if rule violated
logp_violation_matrix <- matrix(unlist(logp_violation)) #making amtrix from list to cbind later
mw_violation <- lapply(mw, function(x) {if(x>500) return(1) else return(0)})#figure if rule violated
mw_violation_matrix <- matrix(unlist(mw_violation)) #making amtrix from list to cbind later
violations <- cbind(smiles_matrix, hba1_violation_matrix, hbd_violation_matrix, logp_violation_matrix, mw_violation_matrix)
colnames(violations) <- c("SMILE", "Hydrogen Bond Donors > 5", "Hydrogen Bond Acceptors > 10", "Octanol/Water Partition Coefficient (logP) > 5", "Molecular Weight > 500")
print(violations)
print(stats)
write.table(violations, "druglikeness-violations.csv", row.names=FALSE, col.names=TRUE, sep=",")
write.table(stats, "druglikeness-stats.csv", row.names=FALSE, col.names=TRUE, sep=",")
