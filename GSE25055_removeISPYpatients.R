#!/usr/bin/env Rscript

#require(utils)
#print(R.Version())
#####modified ArrayExpress example from Purvesh

load("GSE25055_GPL96.RData")

ISPY = read.delim(header = TRUE, "GSE25055_ISPY patients.txt")
ISPY_GSMID = ISPY$GSMID
  ISPYcount =1; count = 1; 
  for (i in 1:310){
  if(GSE25055_GPL96$GSMID[i] != as.character(ISPY_GSMID[ISPYcount])){ 
    expr_corrected[,count] = GSE25055_GPL96$expr[,i]; count = count+1
    }
  else{ISPYcount=ISPYcount+1}
  }
