#!/usr/bin/env Rscript

#require(utils)
#print(R.Version())
#####modified ArrayExpress example from Purvesh


#clean up data space
rm(list=ls())
gc()

require("Affy")
require(gcrma)

array = "gpl_96"

#removed files that were I-SPY patients (already in GSE22226)

filePath = "~/expr_breastCancerDatasets/GSE25055/patientFiles"
fileNames = list.files(filePath)
print(length(fileNames))
print(fileNames)

#can't specify a path with ReadAffy so change directory
setwd("~/expr_breastCancerDatasets/GSE25055/patientFiles")

rawExpressionSet = ReadAffy(filenames=fileNames)
names(rawExpressionSet)
setwd("~/expr_breastCancerDatasets/GSE25055")

probes = read.delim(file = "~/GeneSymbol/GSE25065_GPL96_probes.txt", header = TRUE)
print(colnames(probes))
#probeSet = probes[,1]
#GeneSymbol = probes[,2]

probeSet = featureNames(rawExpressionSet)

Pheno = read.delim(file= "GSE25055_clinical_noISPY.txt", header = TRUE)

#background correct w/ rma, norm and log2 with quantiles using gcrma
dataNorm = gcrma(rawExpressionSet)

source(file = "~/global_R_code/db_functions.r")
con = connect.db("ywrfc09", "aveelyau05", host="buttelab-db1", dbname="annot_gpl")
dataMapping = mapProbesToGenes.db(con, probeSet, map.to="Symbol", 
gpl.table=array)
dbDisconnect(con)

GeneSymbol = dataMapping[probeSet, "gene"]

#clean up GSMID
require("limma")
GSM = strsplit2(fileNames,"_")
GSM = GSM[,1]
GSM = substring(GSM, 4)
 

colnames(exprs(dataNorm)) = GSM

expr = exprs(dataNorm)

GSE25055_GPL96 = list(expr=expr, Pheno=Pheno, GSMID = GSM,ReporterID = 
probeSet,GeneSymbol=GeneSymbol)
boxplot(GSE25055_GPL96$expr, las=3, main="GSE25055_GPL96")
 
#manual check on length of keys, expr needs to have this many rows
print(length(GSE25055_GPL96$GeneSymbol))
#manual check on number of patients I have, expr needs to have this many cols
print(length(GSE25055_GPL96$GSMID))
#manual check that expr is data frame with genes as rows and sample values as cols 
print(dim(GSE25055_GPL96$expr))


#check GSMID
print(GSE25055_GPL96$GSMID[1:10])

#check GeneSymbol
print(GSE25055_GPL96$GeneSymbol[300:330])

save(GSE25055_GPL96, file="GSE25055_GPL96.RData")
