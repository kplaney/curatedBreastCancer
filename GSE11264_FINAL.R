#!/usr/bin/env Rscript

#require(utils)
#print(R.Version())
#####modified ArrayExpress example from Purvesh
library("RMySQL")

#look up array first
array = "gpl_4819"
#read in the phenotype (i.e. clinical information) separately from the .srdf file
#(every Array Express experiment should have a separate srdf file)
#you need to have a file with just the probes in a column in a separate directory.
Pheno =  read.delim("E-GEOD-11264.sdrf.txt", sep="\t", header=TRUE, dec = ".", quote="\"")

#weird...below command not working!
#Probes_full = read.delim("GPL4819_A-GEOD-4819.adf_clean.txt", sep = "t\", header = TRUE)
#print(colnames(Probes_full))
#print(Probes_full[1,])
#print(Pheno[1,])
nameCol = 15
probes =  read.delim("Probes.txt", header = TRUE)
#dim of probes below should match numProbes pulled from an example patient file!
numProbes = dim(probes)[1]
#look at all patient files, and extract the GSMIDs:

fileNames =  list.files("~/expr_breastCancerDatasets/GSE11264/patientFiles")
#check what the file names look like
print(fileNames[1])
numFiles = length(fileNames)
GSM = strsplit(as.character(Pheno[,nameCol]),".txt")
print(GSM[1])
GSM = strsplit(as.character(GSM), ".CEL")
GSM = strsplit(as.character(GSM), ".gpr")
print(GSM[1])
GSM = strsplit(as.character(GSM),"_")
print(GSM[1])
#we only want the first word, the GSMID, in each row
cutoff_letter = "s"
print(length(GSM[1]))
GSM = substring(GSM,4)
print(GSM)
print(length(GSM))
print(numFiles)
print(numProbes)
print(length(probes[,1]))


Expr =  matrix(data = NA, nrow=numProbes, ncol=numFiles, dimnames = 
list(probes[,1],GSM));


source("~/global_R_code/db_functions.r")
con = connect.db(username ="ywrfc09", password = "aveelyau05",host="buttelab-db1", dbname="annot_gpl")
dataMapping = mapProbesToGenes.db(con, probes[,1], map.to="Symbol", gpl.table=array)
dbDisconnect(con)
print(attributes(dataMapping))

#make sure there aren't extraneous rows in the files before importing.
#must know ahead of time which column you want to pull out!
#ALSO: check ahead of time that your probe numbers match the order of what you're pulling out of each patient file!
#column we want in each patient file - ususally lowess normalized value.
#if don't know, open up the .srdf (easy to read in excel) and find out the column value
dataCol= 2;
setwd("~/expr_breastCancerDatasets/GSE11264/patientFiles")
for (i in 1:numFiles) {
  patientFile <- read.delim(fileNames[i], sep="\t", header=TRUE, dec = ".", quote="\"")  
  Expr[,i] <-patientFile[,dataCol];
}

#when I did a boxplot, had negative values, and didn't look really normalized.  was log2 though
#so do it yourself.

quantileNormalize <- function(dataMatrix){

  require(preprocessCore)
  tempData = normalize.quantiles(dataMatrix)
  rownames(tempData) = rownames(dataMatrix)
  colnames(tempData) = colnames(dataMatrix)
  
  #adjust code so that it's all greater than 0
  if(min(tempData,na.rm=T) <0){
    tempData = tempData - min(tempData, na.rm=T)
  } else{
    tempData = tempData + min(tempData, na.rm=T)

 }
  return(tempData)

}

expr = quantileNormalize(Expr)

#then log2( ) data if not already! here it is.

setwd("~/expr_breastCancerDatasets/GSE11264")
GSE11264_GPL4819 <- list(expr = expr, GSMID = GSM, GeneSymbol = dataMapping[probes[,1],"gene"]);
print(attributes(Pheno))
print(GSE11264_GPL4819$GSMID)
print(length(GSE11264_GPL4819$GeneSymbol))
print(dim(GSE11264_GPL4819$expr))
print(numFiles)
print(numProbes)
save(GSE11264_GPL4819, probes, Pheno, file = "GSE11264_GPL4819.RData")
boxplot(GSE11264_GPL4819$expr, las = 3, main = "GSE11264_GPL4819")



