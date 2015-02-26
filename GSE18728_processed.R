#!/usr/bin/env Rscript

#require(utils)
#print(R.Version())
#####modified ArrayExpress example from Purvesh

require("ArrayExpress")
require(gcrma)

array = "gpl_570"
procData = getAE("E-GEOD-18278", type="processed") #download all files from 
attributes(procData)
cnames = getcolproc(procData)
print(cnames)

q()
#rawExpressionSet = magetab2bioc(files = rawData, rawcol = "GenAu-MB-BCELL-v4.0")
#Pheno = pData(phenoData(rawExpressionSet)) #get all the phenodata and save for later, even though don't need it now

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
colnames(exprs(dataNorm)) = GSM
GSE18728_GPL570 = list(expr=exprs(dataNorm), GSMID = GSM,GeneSymbol=dataMapping[rownames(exprs(dataNorm)),"gene"])
boxplot(dataNorm$expr, las=3, main="GSE18728_GPL570") 
#manual check on length of keys, expr needs to have this many rows
length(GSE18728_GPL570$keys)
#manual check on number of patients I have, expr needs to have this many cols
length(GSE18728_GPL570$GSMID)
#manual check that expr is data frame with genes as rows and sample values as cols 
dim(GSE18728_GPL570$expr)

print(GSE187278_GPL570$GSMID)

#look at pheno - which column is GSMIDs?
attributes(Pheno)

print(GSE18728_GPL570$GeneSymbol[1:20])
save(GSE187278_GPL570, Pheno, file="GSE18728_GPL570.RData")
