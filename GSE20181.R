#!/usr/bin/env Rscript

#require(utils)
#print(R.Version())
#####modified ArrayExpress example from Purvesh


#clean up data space
gc()

require("ArrayExpress")
require(gcrma)

array = "gpl_96"
rawData = getAE("E-GEOD-32072", type="raw") #download all files from array 

#magetab2bioc takes log2 for you
rawExpressionSet = magetab2bioc(files = rawData)
Pheno = pData(phenoData(rawExpressionSet)) #get all the phenodata and save for later, even though don't need it now

#background correct w/ rma, norm with quantiles using gcrma
dataNorm = gcrma(rawExpressionSet)
source(file = "~/global_R_code/db_functions.r")
con = connect.db("ywrfc09", "aveelyau05", host="buttelab-db1", dbname="annot_gpl")
dataMapping = mapProbesToGenes.db(con, rownames(exprs(dataNorm)), map.to="Symbol", gpl.table=array)
dbDisconnect(con)
require("limma")
GSM = strsplit2(colnames(exprs(dataNorm)),".CEL")
GSM = GSM[,1]
GSM = substring(GSM,4) 
colnames(exprs(dataNorm)) = GSM
GSE32072_GPL96 = list(expr=exprs(dataNorm), Pheno = Pheno, GSMID = GSM,ReporterID = 
rownames(exprs(dataNorm)),GeneSymbol=dataMapping[rownames(exprs(dataNorm)),"gene"])
boxplot(GSE32072_GPL96$expr, las=3, main="GSE320732_GPL96") 
#manual check on length of keys, expr needs to have this many rows
length(GSE32072_GPL96$GeneSymbol)
#manual check on number of patients I have, expr needs to have this many cols
length(GSE32072_GPL96$GSMID)
#manual check that expr is data frame with genes as rows and sample values as cols 
dim(GSE32072_GPL96$expr)

print(GSE32072_GPL96$GSMID[1:10])

#look at pheno - which column is GSMIDs?
attributes(Pheno)

print(GSE32072_GPL96$GeneSymbol[1:20])

save(GSE32072_GPL96, file="GSE32072_GPL96.RData")
