#the new Bioconductor checks require protocolData to be added:

data("curatedBreastDataExprSetList")
curatedBreastDataExprSetList  <- lapply(curatedBreastDataExprSetList, function(x){
labelDescription <- rep("Breast cancer human tumor tissue sample", nrow(pData(x)))
labelDescription <- data.frame(labelDescription) 
rownames(labelDescription) <- rownames(pData(x))
protocolData(x) <- AnnotatedDataFrame(labelDescription)
return(x)
})

table(sapply(curatedBreastDataExprSetList, function(x)
 isTRUE(validObject(x, test=TRUE))))

#test - does processing still work?
proc_curatedBreastDataExprSetList <- processExpressionSetList(
 exprSetList=curatedBreastDataExprSetList, 
 outputFileDirectory = "~/Desktop/", numTopVarGenes=5000)

#NOTE: need to provide full file path.
#save(curatedBreastDataExprSetList,file = "curatedBreastDataExprSetList.rda",compress='xz')
