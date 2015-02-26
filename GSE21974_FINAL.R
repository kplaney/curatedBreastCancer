#!/usr/bin/env Rscript

#require(utils)
#print(R.Version())
#####modified ArrayExpress example from Purvesh

#clear workspace
rm(list=ls())

#script for processed expression data from GEO
#this data is already quantile normalized

#what array the study used
array = "gpl_6480"
#download my cleaned-up series matrix from GEO. first column is probeID
data = read.delim(header = TRUE, file = "GSE21974_preTreat_expression.txt")
expr = data[,2:dim(data)[2]]
print(which(is.na(expr)))

ProbeID = data[,1]

GSMID = substring(colnames(expr), 4)
print(GSMID[1:10])

#looked-already all nonzero values

print(dim(expr))
summary(expr)


#load Purvesh's database functions db_functions.r
source(file = "~/global_R_code/db_functions.r")
con = connect.db("ywrfc09", "aveelyau05", host="buttelab-db1", dbname="annot_gpl")

#extract genes from your probe list
dataMapping = mapProbesToGenes.db(con, ProbeID, map.to="Symbol", gpl.table=array)
dbDisconnect(con)

summary(dataMapping)
print((dataMapping[ProbeID,"gene"])[1:20])

reporterID = ProbeID

geneSymbol = dataMapping[ProbeID, "gene"]
 
boxplot(expr, las = 3, main = "GSE21974_GPL6480_pre-treatment") 

setwd("~/expr_breastCancerDatasets/GSE21974")
Pheno = read.delim(header = TRUE, file = "GSE21974_preTreat_clinical.txt")

GSE21974_GPL6480= list(expr=expr, GSMID= GSMID, geneSymbol = geneSymbol, 
reporterID = reporterID, Pheno = Pheno)

save(GSE21974_GPL6480, file = "GSE21974_GPL6480_FINAL.RData")

#check dim
print(length(GSE21974_GPL6480$geneSymbol))
print(length(GSE21974_GPL6480$reporterID))
print(length(GSMID))

#row length of expr should equal row length of gene symbol!!
print(dim(GSE21974_GPL6480$expr))                                 
