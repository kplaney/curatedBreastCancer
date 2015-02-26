#!/usr/bin/env Rscript

#require(utils)
#print(R.Version())
#####modified ArrayExpress example from Purvesh
library("RMySQL")

array = "gpl_180"
#read in the phenotype (i.e. clinical information) separately from the .srdf file
#(every Array Express experiment should have a separate srdf file)
#you need to have a file with just the probes in a column in a separate directory.
Pheno =  read.delim("E-GEOD-61.sdrf.txt", sep="\t", header=TRUE, dec = ".", quote="\"")
probes_full = read.delim("GPL180_A-GEOD-180_clean.adf.txt", sep = "\t", header = TRUE, dec = ".", quote = "\"")
print(colnames(probes_full))
print(probes_full[1:5,])
print(Pheno[1,])
nameCol = 1
probes =  read.table("GPL180_genbank_probes.txt", header = TRUE);
#dim of probes below should match numProbes pulled from an example patient file!
numProbes = dim(probes)[1]
#look at all patient files, and extract the GSMIDs:

fileNames =  list.files("~/expr_breastCancerDatasets/GSE61/patientFiles");
#check what the file names look like
print(fileNames[1])
numFiles = length(fileNames)
print(Pheno[,nameCol])
GSM = strsplit(as.character(Pheno[,nameCol]),".txt")
GSM = strsplit(as.character(GSM), " ")
GSM = strsplit(as.character(GSM), ".CEL")
GSM = strsplit(as.character(GSM), ".gpr")
GSM = strsplit(as.character(GSM),"_")
GSM = substring(GSM,4)
print(GSM)
print(length(GSM))
print(numFiles)
print(numProbes)
#print(length(probes[,1]))
#look in a sample patient file. which one is the normalized column?
probeCol = 11
print(probes[1,probeCol])
Expr =  matrix(data = NA, nrow=numProbes, ncol=numFiles, dimnames = 
list(probes[,probeCol],GSM));

con = connect.db(username ="ywrfc09", password = "aveelyau05",host="buttelab-db1", dbname="annot_gpl")
dataMapping = mapProbesToGenes.db(con, probes[,1], map.to="Symbol", gpl.table=array)
dbDisconnect(con)
#print(attributes(dataMapping))

#make sure there aren't extraneous rows in the files before importing.
#must know ahead of time which column you want to pull out!
#ALSO: check ahead of time that your probe numbers match the order of what you're pulling out of each patient file!
#column we want in each patient file - ususally lowess normalized value.
#if don't know, open up the .srdf (easy to read in excel) and find out the column value
dataCol= 2;
setwd("~/expr_breastCancerDatasets/GSE61/patientFiles")
for (i in 1:numFiles) {
  patientFile <- read.delim(fileNames[i], sep="\t", header=TRUE, dec = ".", quote="\"")  
  Expr[,i] <-patientFile[,dataCol];
}
setwd("~/expr_breastCancerDatasets/GSE61")
GSE61_GPL180 <- list(expr = Expr, GSMID = GSM, keys = 
dataMapping[probes[,1],"gene"], KPkeys = probes_full);
print(attributes(Pheno))
print(GSE61_GPL61$GSMID)
print(dim(GSE61_GPL61$keys))
print(dim(GSE61_GPL61$expr))
print((GSE61_GPL61$expr)[1:10,1:10])
print(numFiles)
print(numProbes)
save(GSE61_GPL61, probes, probes_full, Pheno, file = "GSE61_GPL180.RData")
#print(GSE61_GPL180$keys)
#print(probes)
#check that data is normalized, and is log2
#check with a box plot on your own computer later! 
#GSE22358boxplot <- boxplot(GSE22358$expr, las=3, main="E-GEOD-22358")




