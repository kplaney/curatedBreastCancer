#!/usr/bin/env Rscript

#require(utils)
#print(R.Version())
#####modified ArrayExpress example from Purvesh
require("limma")
require("RmySQL")

#script for .gpr files

#file path to separate directory of downloaded .gpr files

filePath = 
("~/expr_breastCancerDatasets/GSE21974/GSE21974_RAW")

#what array the study used
array = "gpl_6480"

fileNames = list.files(filePath)

imageAnalysisProgram = "spot"
RG = read.maimages(fileNames, source = imageAnalysisProgram, path = 
filePath, columns =list(Rf="rMedianSignal", Gf="gMedianSignal", Rb 
="rBGMedianSignal", Gb="gBGMedianSignal"))
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
print(GSMID[1:10])


print(rownames(MA.p$A)[1:10])

expr = MA.p$A
summary(expr)

print(names(MA.p$genes))
rownames(expr) = (MA.p$genes)$ID
print(rownames(expr)[1:10])


#load Purvesh's database functions db_functions.r
#source(file = "~/global_R_code/db_functions.r")
#con = connect.db("ywrfc09", "aveelyau05", host="buttelab-db1", dbname="annot_gpl")

#extract genes from your probe list
dataMapping = mapProbesToGenes.db(con, (MA.p$genes)$ID, map.to="Symbol", gpl.table=array)
dbDisconnect(con)

summary(dataMapping)
print((dataMapping[rownames(MA.p$genes)$ID,"gene"])[1:20])
q()

#get Gene Symbol from your database

#setwd("~/GeneSymbol")
#load("GPL6848_GeneSymbol.RData")
#pull out the Gene Symbol list you previously created

#why I can't also use "genes" from limma: sometimes lists multiple genes
#on one line so you should clean up the genes text beforehand
q()
names(annot)

#double-check: are all of the probes actually used in this specific study?

#make sure na.rm below removes NA from BOTH columns...

if(all(rownames(expr) == annot$ID, na.rm = TRUE)){
 reporterID = annot$ID
 geneSymbol = annot$GeneSymbol
 print("your targetIDs and the study targetIDs match up")

}else{

 probeCount =1
 reporterID = array(data = NA, dim = dim(expr)[1])
 geneSymbol = array(data = NA, dim = dim(expr)[1])

  for(j in 1:length(annot$ID)){

   if((annot$ID)[j] == rownames(expr)[probeCount]){

   reporterID[probeCount] = annot$ID[j]
   geneSymbol[probeCount] = annot$GeneSymbol[j]
   probeCount = probeCount +1

   }
 }

}

#now double-check again

print(all(rownames(expr) == reporterID))
#if you get NA, that's fine - just means some rows are NA,
#but you still didn't get any false values.

boxplot(expr, las = 3, main = "GSE18378_GPL6848_pre-treatment") 

setwd("~/expr_breastCancerDatasets/GSE18378")
Pheno = read.delim(header = TRUE, file = 
"GSE18378_GPL6848_preTreat_clinical.txt")

GSE18378_GPL6848 list(expr=expr, GSMID= GSMID, geneSymbol = 
geneSymbol, reporterID = reporterID,  Pheno = Pheno)

save(GSE18378_GPL6848, file = "GSE18378_GPL6848-preTreat.RData")

#check dim
print(length(GSE18378_GPL6848$geneSymbol))
print(length(GSE18378_GPL6848$reporterID))
print(length(GSMID))

#row length of expr should equal row length of gene symbol!!
print(dim(GSE18378_GPL6848$expr))                                 
