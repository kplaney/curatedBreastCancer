#!/usr/local/R/bin/Rscript --vanilla

#require(utils)
#print(R.Version())
##### 

rm(list=ls())
gc()

#example of study that only released processed data, on two arrays

array = "gpl_4133"

#I already downloaded the processed text files via getAE( ) type = "processed"
filePath = "~/expr_breastCancerDatasets/GSE22226/patientFiles_GPL4133"

#one file for each patient
fileNames = list.files(filePath, full.names = TRUE)
shortFileNames = list.files(filePath)

testFile = read.delim(fileNames[1], header = TRUE)

#first column is probe list
probeSet = testFile[,1]


#Purvesh's data mapping didn't work with these probes 
#so created my own list from GEO site
#load("~/GeneSymbol/GPL4133_GeneSymbol.RData")
#names(annot)
#probeSet = annot$ID
#GeneSymbol = annot$GeneSymbol

#check GeneSymbol output with some random rows
#print(GeneSymbol[300:330])

#loop through and pull out expression values from files
expr = matrix(data =NA, nrow = dim(testFile)[1], ncol = length(fileNames))

for (i in 1:length(fileNames)){

        tempFile = read.delim(fileNames[i], header = TRUE)

	#pull out second column with expression values
	expr[,i] = tempFile[,2]

}


#this data is already normalized with MAS 5.0, all and log2.
#make positive? naw keep original values for now.
#expr = expr -min(expr, na.rm = TRUE)

#get GSMIDs
require("limma")
GSMID = strsplit2(shortFileNames, "_")
GSMID = GSMID[,1]
GSMID = substring(as.character(GSMID),4)

#double-check
print(GSMID[1:10])


#already plotted to check
boxplot(expr, las=3, main = "GSE22226_GPL4133")


GSE22226_GPL4133 = list(expr=expr, keys = probeSet, GSMID=GSMID)
save(GSE22226_GPL4133, file = "GSE22226_GPL4133.RData")

print(dim(GSE22226_GPL4133$expr))
print(length(GSE22226_GPL4133$keys))
print(length(GSE22226_GPL4133$GSMID))

