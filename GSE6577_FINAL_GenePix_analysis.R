#!/usr/bin/env Rscript

#require(utils)
#print(R.Version())
#####modified ArrayExpress example from Purvesh
require("limma")
require("RmySQL")

#script for .gpr files

#file path to separate directory of downloaded .gpr files

filePath = 
("~/expr_breastCancerDatasets/GSE6577/GSE6577_RAW")

#what array the study used
array = "gpl_3883"

fileNames = list.files(filePath)

imageAnalysisProgram = "genepix"
RG = read.maimages(fileNames, source = imageAnalysisProgram, path = filePath)
#background correct using suggested method for differential expression values in Limma R guide
#adjust foreground adaptively for background intensities and results in strictly positive
#adjust intenstities (not negative or zero corrected values)
backMethod = "normexp"
#offset: damps the variation of the log ratios for very low intensity spots towards zero
ofSet = 50
RGb = backgroundCorrect(RG, method = backMethod,offset = 50)

#assuming print-tip loss. if agilent array: add method = "loess".
#if small arrays with less than approx 150 spots per print-tip: use global "loess"
#normalization or "robustspline" normalization
#I use "loess" to be on the safe side...quality of some of these studies is questionable
#plus, we're using loess normalization for most of our agilent studies.

MA = normalizeWithinArrays(RGb, method = "loess")

#now normalize between patients. limma defaults to Aquantile for 2 colors (leaves log values unchanged)
#and quantile for one colors (we usually have 2 colors)

#I used the normalized cyclic loess "affy" method, as a lot of our affy arrays use that.

MA.p = normalizeBetweenArrays(MA, method = "cyclicloess", cyclic.method = 
"affy")

GSMID = substring(rownames(MA.p$targets), 4)
print(GSMID[1])


probes = (MA.p$genes)$ID

#add row number to make sure we align everything
genes = list(genes = (MA.p$genes)$Name, rowIndex = (MA.p$genes)$Row)

expr = MA.p$M

#load Purvesh's database functions db_functions.r
source(file = "~/global_R_code/db_functions.r")
con = connect.db("ywrfc09", "aveelyau05", host="buttelab-db1", dbname="annot_gpl")

#extract genes from your probe list
dataMapping = mapProbesToGenes.db(con, probes, map.to="Symbol", gpl.table=array)
dbDisconnect(con)

summary(dataMapping)

GSE6577_GPL3883= list(expr=expr, GSMID= GSMID,keys=dataMapping, geneList = genes, probes = probes)

save(GSE6577_GPL3883, file = "GSE6577_GPL3883.RData")

#check dim
print(length((GSE6577_GPL3883$genes)
print(length(probes))
print(length(GSMID))
print(dim(GSE6577_GPL3883$expr))
              
#box plot to check your work
boxplot(GSE6577_GPL3883$expr, las=3, main="E-MTAB-6008") #run on local machine
                                  

