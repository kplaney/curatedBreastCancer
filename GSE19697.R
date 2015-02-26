#!/usr/bin/env Rscript

#require(utils)
#print(R.Version())
#####modified ArrayExpress example from Purvesh

rm(list=ls())
gc()

require("affy")
require("gcrma")

#script for non AE, affy (.CEL) files:use Affy package

#file path to separate directory of downloaded .CEL files

filePath = 
("~/expr_breastCancerDatasets/GSE19697/GSE19697_RAW")

#what array the study used
array = "gpl_570"

fileNames = list.files(filePath, full.names=TRUE)
fileNamesShort = list.files(filePath)

setwd("~/expr_breastCancerDatasets/GSE19697/GSE19697_RAW")
data = ReadAffy(filenames = fileNamesShort)

dataNorm = gcrma(data)

expr = exprs(dataNorm)
setwd("~/expr_breastCancerDatasets/GSE19697")


data = read.delim("GSE19697_expression.txt", header = TRUE)
#arrays = data[,2:(dim(data)[2])]
probeSet = data[,1]


print(summary(expr))
print(dim(expr))



require("limma")
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
dataMapping = mapProbesToGenes.db(con, probeSet, map.to="Symbol", gpl.table=array)
dbDisconnect(con)

summary(dataMapping)
print((dataMapping[probeSet,"gene"])[1:30])

#get Gene Symbol from your database

#setwd("~/GeneSymbol")

#I created a special Gene Symbol list for this datasets - the expr matrix 
#included all of the NA values. bummer!

#load("GPL6480_GSE21794special_GeneSymbol.RData")
#pull out the Gene Symbol list you previously created

#why I can't also use "genes" from limma: sometimes lists multiple genes
#on one line so you should clean up the genes text beforehand

#names(annot)

#double-check: are all of the probes actually used in this specific study?

#make sure na.rm below removes NA from BOTH columns...
#reporterID = annot$ProbeName
#geneSymbol = annot$GeneSymbol

#if(all(rownames(expr) == annot$ID, na.rm = TRUE)){
# reporterID = annot$ID
 #geneSymbol = annot$GeneSymbol
 #print("your targetIDs and the study targetIDs match up")

#}else{

 #probeCount =1
 #reporterID = array(data = NA, dim = dim(expr)[1])
 #geneSymbol = array(data = NA, dim = dim(expr)[1])

  #for(j in 1:length(annot$ID)){

   #if((annot$ID)[j] == rownames(expr)[probeCount]){

   #reporterID[probeCount] = annot$ID[j]
   #geneSymbol[probeCount] = annot$GeneSymbol[j]
   #probeCount = probeCount +1

   #}
 #}

#}

#now double-check again

#print(all(rownames(expr) == reporterID))
#if you get NA, that's fine - just means some rows are NA,
#but you still didn't get any false values.
#However make sure it's not ALL NAs!

boxplot(expr, las = 3, main = "GSE19697_GPL570") 

setwd("~/expr_breastCancerDatasets/GSE19697")
Pheno = read.delim(header = TRUE, file = "GSE19697_clinical.txt")

GSE19697_GPL570= list(expr=expr, GSMID= GSMID, geneSymbol = 
dataMapping[probeSet,"gene"], reporterID = probeSet, Pheno = Pheno)

save(GSE19697_GPL570, file = "GSE19697_GPL570.RData")

#check dim
print(length(GSE19697_GPL570$geneSymbol))
print(length(GSE19697_GPL570$reporterID))
print(length(GSMID))

#row length of expr should equal row length of gene symbol!!
print(dim(GSE19697_GPL570$expr))                                 

