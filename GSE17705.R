#!/usr/bin/env Rscript

#require(utils)
#print(R.Version())
#####modified ArrayExpress example from Purvesh

rm(list=ls())
gc()

#standard Array Express script

require("ArrayExpress")
require(gcrma)

array = "gpl_96"
rawData = getAE("E-GEOD-17705", type="raw") #download all files from 
rawExpressionSet = magetab2bioc(files = rawData)
Pheno = pData(phenoData(rawExpressionSet)) #get all the phenodata and save for later, even though don't need it now

dataNorm = gcrma(rawExpressionSet)


#get Gene Symbols using Purvesh's code
source("~/global_R_code/db_functions.r")
con = connect.db("ywrfc09", "aveelyau05", host="buttelab-db1", dbname="annot_gpl")
dataMapping = mapProbesToGenes.db(con, rownames(exprs(dataNorm)), map.to="Symbol", gpl.table=array)
dbDisconnect(con)

require("limma")
GSM = strsplit2(colnames(exprs(dataNorm)),".CEL")
GSM = GSM[,1]
GSM = strsplit2(GSM, "_")
GSM = GSM[,1]
GSM = substring(GSM,4) 

#check GSM
print(GSM[1:10])

colnames(exprs(dataNorm)) = GSM

GSE17705_GPL96 = list(expr=exprs(dataNorm), GSMID = 
GSM,GeneSymbol=dataMapping[rownames(exprs(dataNorm)),"gene"])


boxplot(GSE17705_GPL96$expr, las=3, main="GSE17705_GPL96") 
#manual check on length of keys, expr needs to have this many rows
length(GSE17705_GPL96$GeneSymbol)
#manual check on number of patients I have, expr needs to have this many cols
length(GSE17705_GPL96$GSMID)
#manual check that expr is data frame with genes as rows and sample values as cols 
dim(GSE17705_GPL96$expr)

print(GSE17705_GPL96$GSMID)

#look at pheno - which column is GSMIDs?
attributes(Pheno)

print(GSE17705_GPL96$GeneSymbol[1:20])
save(GSE17705_GPL96, Pheno, file="GSE17705_GPL96.RData")
