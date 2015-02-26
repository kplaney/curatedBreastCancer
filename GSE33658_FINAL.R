#!/usr/bin/env Rscript

#require(utils)
#print(R.Version())
#####modified ArrayExpress example from Purvesh

#Affy/.CEL files not in AE

#clean workspace
rm(list=ls())
require("affy")
require("gcrma")
#script for non AE, affy (.CEL) files.

#file path to separate directory of downloaded .CEL files

filePath = 
("~/expr_breastCancerDatasets/GSE33658/patientFiles")

#what array the study used
array = "gpl_570"

fileNames = list.files(filePath, full.names=TRUE)

#we have an Affymetrix array that isn't in Array Express. i.e .CEL files
#so use the affy package
require("affy")

#odd but just setting the file path in ReadAffy isn't working. do it before.
("~/expr_breastCancerDatasets/GSE33658/patientFiles")
print(fileNames)
arrays = ReadAffy(filenames = fileNames)
features = featureNames(arrays)
print(features[1:10])
summary(arrays$featureData)
probeSet = features

#background correct and normalize using gcrma( ) like Purvesh does with AE datasets
#this uses quantile normalization
normArrays = gcrma(arrays, normalize = TRUE)

#gcrma outputs an ExpressionSet. a class-specific function
#use exprs( ) to get a summary expression value for each probe/row
expr = exprs(normArrays)
print(summary(expr))
print(dim(expr))


#get your unique GEO GSM IDs
require("limma")
fileNamesShort = list.files(path=filePath)
GSMID = substring(fileNamesShort, 4)
GSMID = strsplit2(as.character(GSMID), ".CEL")
print(attributes(GSMID))
print(names(GSMID))
GSMID = GSMID[,1]
print(GSMID[1:10])


#load Purvesh's database functions db_functions.r
source(file = "~/global_R_code/db_functions.r")
con = connect.db("ywrfc09", "aveelyau05", host="buttelab-db1", dbname="annot_gpl")

#extract genes from your probe list
dataMapping = mapProbesToGenes.db(con, probeSet, map.to="Symbol", 
gpl.table=array)
dbDisconnect(con)

summary(dataMapping)
print((dataMapping[probeSet,"gene"])[1:30])

boxplot(expr, las = 3, main = "GSE33658_GPL570_pre-treatment") 

setwd("~/expr_breastCancerDatasets/GSE33658")
Pheno = read.delim(header = TRUE, file = "GSE33658_preTreat_clinical.txt")

GSE33658_GPL570= list(expr=expr, GSMID= GSMID, ReporterID = 
probeSet, GeneSymbol = 
dataMapping[probeSet,"gene"], reporterID = probeSet, Pheno = Pheno)

save(GSE33658_GPL570, file = "GSE33658_GPL570-preTreat_FINAL.RData")

#check dim
print(length(GSE33658_GPL570$geneSymbol))
print(length(GSE33658_GPL570$reporterID))
print(length(GSMID))

#row length of expr should equal row length of gene symbol!!
print(dim(GSE33658_GPL570$expr))                                 

