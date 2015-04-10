####

createS4exprSet <- function(expr,phenoData,featureData){
  
  #search for class.
  if((!missing(phenoData))&&(!missing(featureData))){
    
    
    exprSet <- new("ExpressionSet", assayData = assayDataNew(exprs=
                   new("matrix")), phenoData = new("AnnotatedDataFrame"), 
                   featureData = new("AnnotatedDataFrame"), 
                   experimentData = new("MIAME"), annotation = character(0))
    
    exprSet@assayData <- assayDataNew(exprs = expr)
    
    rownames(featureData) <- rownames(expr)
    featureData <-  new("AnnotatedDataFrame", data=as.data.frame(featureData))
    exprSet@featureData <- featureData
    
    rownames(phenoData) <- colnames(expr)
    
    phenoData <-  new("AnnotatedDataFrame", data=as.data.frame(phenoData))
    exprSet@phenoData <- phenoData
    #you can everything in your "new" implementation of exprSet, or just add 
    #the expression data (must add these this in new() 
    #or else can't add it later!)
    #sample names must be the same
    #how you would access the expression data: exprSet@assayData$expr
    #add the other data
    
    
  }else if((missing(phenoData))&&(!missing(featureData))){
    
    rownames(featureData) <- rownames(expr)
    featureData <-  new("AnnotatedDataFrame", data=as.data.frame(featureData))
    
    exprSet <- new("ExpressionSet", assayData = 
                     assayDataNew(exprs=new("matrix")), 
                   phenoData = new("AnnotatedDataFrame"), 
                   featureData = new("AnnotatedDataFrame"), 
                   experimentData = new("MIAME"), annotation = character(0))
    
    exprSet@assayData <- assayDataNew(exprs = expr)
    
    exprSet@featureData <- featureData
    
  }else if ((!missing(phenoData))&&(missing(featureData))){
    
    rownames(phenoData) <- colnames(expr)
    
    phenoData <-  new("AnnotatedDataFrame", data=as.data.frame(phenoData))
    
    exprSet <- new("ExpressionSet", 
                   assayData = assayDataNew(exprs=new("matrix")), 
                   phenoData = new("AnnotatedDataFrame"), 
                   featureData = new("AnnotatedDataFrame"), 
                   experimentData = new("MIAME"), annotation = character(0))
    
    exprSet@assayData <- assayDataNew(exprs = expr)
    
    exprSet@phenoData <- phenoData
    
  }else{
    
    exprSet <- new("ExpressionSet", 
                   assayData = assayDataNew(exprs=new("matrix")), 
                   phenoData = new("AnnotatedDataFrame"), 
                   featureData = new("AnnotatedDataFrame"), 
                   experimentData = new("MIAME"), annotation = character(0))
    
    exprSet@assayData <- assayDataNew(exprs = expr)
    
  }
  
  message("ExpressionSet successfully made. 
          Dimensions of expr in AssayData is ",dim(exprSet@assayData$exprs))
  
  #return expression set
  exprSet
  
}

#####
#if have an expressionmatrix list (may only have one element).
#can already have a phenoData name index in each list index,
#or a masterPheno data frame that's for all samples
createExpressionSetList <- function(exprMatrixList,masterPhenoData,
                                    patientKey="GEO_GSMID"){
  
  origNames <- names(exprMatrixList);
  if(length(exprMatrixList)==0){
    
    stop("It looks like your matrix list has no elements!")
    
  }
  
  ExpressionSetList <- list()
  
  for(d in 1:length(exprMatrixList)){
    
    
    if(!is.matrix(exprMatrixList[[d]]$expr)){
      
      warning("This code assumes your expression data in the expr 
              slot is a matrix, as this allows for duplicated gene row names.")
      
    }
    #make an expression object for each exprMatrixList
    
    #add study Name,unique ID for phenoData
    phenoData <- rep(names(exprMatrixList)[d],ncol(exprMatrixList[[d]]$expr))
    
    if(!all(is.null(exprMatrixList[[d]]$phenoData))){
      
      pheno <- exprMatrixList[[d]]$phenoData
      
      #row # is number of patients here.
      phenoData <- cbind(phenoData,pheno)
      
      #user provided a master data list instead - pull out the relevant patients
      #for this loop.
    }else if(!missing(masterPhenoData)){
      
      #usually numeric vs. character not important here
      #(if string is actually all numbers.)
      pData <- masterPhenoData[na.omit(match(colnames(exprMatrixList[[d]]$expr),
                                             masterPhenoData[,patientKey])), ]
      #keeps converting to a list!!! have do add in any potential 
      #NA patients in a roundabout way..
      if(any(is.na(match(colnames(exprMatrixList[[d]]$expr),
                         masterPhenoData[,patientKey])))){
        
        stop(paste0("must have pheno data (even if NA) for each patient.\n",
                    "Error occured for dataset ",d,"\n"))
      }
      
      phenoData <- data.frame(phenoData,pData)
      
    }
    
    #Row names of phenoData must match column names of the 
    #matrix / matricies in expr.
    rownames(phenoData) <- colnames(exprMatrixList[[d]]$expr)
    #add a variable for study name.
    colnames(phenoData)[1] <- c("datasetName")
    
    #store row names in feature data - see below how will need to update 
    #row names of assay data to store featureData.
    featureData <- rownames(exprMatrixList[[d]]$expr)
    
    
    if(!all(is.null(exprMatrixList[[d]]$featureData))){
      
      featureData <- data.frame(featureData,exprMatrixList[[d]]$featureData)
      
    }
    
    colnames(featureData)[1] <- "origAssayDataRowName"
    
    #Row names of featureData data frane must match row names of the 
    #matrix / matricies in expr.
    #unfortunately, can't have duplicated rownames in a dataframe...
    #so just make unique IDs for assay data matrix.
    #even probes may have NAs, which will be interpreted 
    #as a duplicated row name.
    rownames(exprMatrixList[[d]]$expr) <- c(1:nrow(exprMatrixList[[d]]$expr))
    rownames(featureData) <- rownames(exprMatrixList[[d]]$expr)
    
    expressionSetList[[d]] <- createS4exprSet(expr=
                              exprMatrixList[[d]]$expr,
                              phenoData=phenoData,featureData=featureData)
    
    
  }
  
  if(length(origNames)==length(expressionSetList)){
    #add in original study names if not NULL
    names(expressionSetList) <- origNames;
    
  }
  
  return(expressionSetList)
  
  }
