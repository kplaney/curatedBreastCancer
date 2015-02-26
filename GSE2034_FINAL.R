#!/usr/bin/env Rscript

#require(utils)
#print(R.Version())
#####modified ArrayExpress example from Purvesh

require("ArrayExpress")
require(gcrma)

array = "gpl_96"
rawData = getAE("E-GEOD-2034", type="raw") #download all files from array express for this dataset
rawExpressionSet = magetab2bioc(files = rawData)
Pheno = pData(phenoData(rawExpressionSet)) #get all the phenodata and save for later, even though don't need it now

dataNorm = gcrma(rawExpressionSet)
expr = exprs(dataNorm)
#plot your data
boxplot(exprs(dataNorm), las= 3, main = "GSE2034_GPL96")
con = connect.db("ywrfc09", "aveelyau05", host="buttelab-db1", 
dbname="annot_gpl")
dataMapping = mapProbesToGenes.db(con, rownames(expr), map.to="Symbol", gpl.table=array)
dbDisconnect(con)
require("limma")
GSM = strsplit2(colnames(exprs(dataNorm)),".CEL")
GSM = GSM[,1]
#GSM = strsplit(GSM, ".CEL")
GSM = substring(GSM,4) 
colnames(exprs(dataNorm)) = GSM
GSE2034_GPL96 = list(expr=expr, GSMID = 
GSM,GeneSymbol=dataMapping[rownames(expr),"gene"], ReporterID = 
rownames(exprs(dataNorm)), Pheno = Pheno)

#manual check on length of keys, expr needs to have this many rows
length(GSE2034_GPL96$GeneSymbol)
#manual check on number of patients I have, expr needs to have this many cols
length(GSE2034_GPL96$GSMID)
#manual check that expr is data frame with genes as rows and sample values as cols 
dim(GSE2034_GPL96$expr)

#look at pheno - which column is GSMIDs?
attributes(Pheno)
print(GSE2034_GPL96$GSMID)
save(GSE2034_GPL96, file = "GSE2034_GPL96.RData")
