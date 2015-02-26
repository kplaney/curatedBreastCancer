#!/usr/bin/env Rscript

#require(utils)
#print(R.Version())
#####modified ArrayExpress example from Purvesh

array = "gpl_6480"

setwd("~/expr_breastCancerDatasets/GSE12665")
procData = 
list.files("~/expr_breastCancerDatasets/GSE12665/patientFiles")
print(procData[1:10])
require("limma")
GSMID = strsplit2(procData, "_")
GSMID = GSMID[,1]
GSMID = substring(as.character(GSMID), 4)
print(GSMID[4])

setwd("~/expr_breastCancerDatasets/GSE12665/patientFiles")
samplePatient = read.delim(procData[1], header = TRUE)
print(dim(samplePatient))
print(colnames(samplePatient))

#load your probe IDs

annot = read.delim(header=TRUE,"~/GeneSymbol/GSE12665_GPL6480_specialAnnotation.txt") 
ReporterID = annot$ID

#GeneSymbol = annot$GeneSymbol

#probes = read.delim("~/GeneSymbol/GPL5325_GeneSymbol.txt", header = TRUE)
expr = matrix(ncol = length(procData), nrow = dim(samplePatient)[1], data = NA, dimnames = list(ReporterID, GSMID))

for(i in 1:length(procData)){
  #skip header/first line - don't want to tack on
  patFile = read.delim(header=FALSE, skip =1, file = procData[i])
  #pull out column you want
  expr[,i] = patFile[,2]
}

setwd("~/expr_breastCancerDatasets/GSE12665")

#get Gene Symbol from Purvesh's dataset

source("~/global_R_code/db_functions.r")

con = connect.db("ywrfc09", "aveelyau05", host = "buttelab-db1", dbname = "annot_gpl")
dataMapping = mapProbesToGenes.db(con, ReporterID, map.to = "Symbol", gpl.table = array)
dbDisconnect(con)

#check do Gene Symbols look OK?
summary(dataMapping)
print((dataMapping[ReporterID, "gene"])[1:30])

GeneSymbol = dataMapping[ReporterID, "gene"]

Pheno = read.delim("E-GEOD-12665.sdrf.txt", sep="\t", header=TRUE, dec = ".", quote="\"")

GSE12665_GPL6480 = list(GeneSymbol = GeneSymbol, ReporterID = 
ReporterID, GSMID = GSMID, expr = expr, Pheno = Pheno)

#check log2???
boxplot(GSE12665_GPL6480$expr, las=3, main = "GSE12665_GPL6480")
 
save(GSE12665_GPL6480, file = "GSE12665_GPL6480.RData")
print(dim(GSE12665_GPL6480$expr))
print(length(GSE12665_GPL6480$GSMID))
print(length(GSE12665_GPL6480$GeneSymbol))
print(length(GSE12665_GPL6480$ReporterID))





