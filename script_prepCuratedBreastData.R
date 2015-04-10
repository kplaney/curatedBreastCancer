


#prep patient data
setwd("~/Dropbox/github/curatedBreastCancer/curatedBreastCancer/")
load("patientData.RData.gzip")

#last few variables are more experimental/processed - remove these.
masterPhenoData <- patientDataFull[,c(1:112)];

dbUniquePatientIDs <- masterPhenoData$dbUniquePatientID;

library('RMySQL')
m <- dbDriver("MySQL")
con<- dbConnect(MySQL(), user="root", dbname="test")
res <- dbSendQuery(conn, "SELECT therapy_chemical_name FROM breastCancer_treatment_dictionary")
drugs<- fetch(res, n=-1)
dbDisconnect(conn)
drugs <- as.vector(t(drugs))

con<- dbConnect(MySQL(), user="root", dbname="test")
res <- dbSendQuery(con, "SELECT * FROM breastCancer_treatment_dictionary")
classes<- fetch(res, n=-1)
dbDisconnect(con)

#COME BACK: add site_prefix here???
#main clincal table: patient_info
treatmentQuery <- paste("select patient_info.dbUniquePatientID as dbUniquePatientID, Y.* FROM patient_info INNER JOIN 
                        
                        ( SELECT  breastCancer_treatment_dictionary.*,Z.study_specific_protocol_number,  
                        Z.study_specific_trt_number, Z.study_ID,Z.neoadjuvant_or_adjuvant, Z.dose, Z.units, Z.time_min, Z.days_cycle1_started_from_day1,
                        Z.day_of_cycle_started,Z.cycle_length_days, Z.num_daysPerCycle, Z.num_timesPerDay,Z.sequential_days,Z.num_cycles FROM breastCancer_treatment_dictionary INNER JOIN
                        
                        ( SELECT X.study_specific_protocol_number,  breastCancer_treatment_ind_therapies.* FROM breastCancer_treatment_ind_therapies 
                        
                        INNER JOIN (SELECT breastCancer_treatment_protocols.study_specific_protocol_number,
                        breastCancer_treatment_protocols.study_ID, breastCancer_treatment_protocols.study_specific_trt_num1, 
                        breastCancer_treatment_protocols.study_specific_trt_num2,  breastCancer_treatment_protocols.study_specific_trt_num3,
                        breastCancer_treatment_protocols.study_specific_trt_num4,breastCancer_treatment_protocols.study_specific_trt_num5,
                        breastCancer_treatment_protocols.study_specific_trt_num6 FROM breastCancer_treatment_protocols ) as X 
                        
                        ON 
                        ( (breastCancer_treatment_ind_therapies.study_ID= X.study_ID) AND 
                        ( (breastCancer_treatment_ind_therapies.study_specific_trt_number=X.study_specific_trt_num1) OR 
                        (breastCancer_treatment_ind_therapies.study_specific_trt_number=X.study_specific_trt_num2) OR
                        (breastCancer_treatment_ind_therapies.study_specific_trt_number=X.study_specific_trt_num3) OR
                        (breastCancer_treatment_ind_therapies.study_specific_trt_number=X.study_specific_trt_num4) OR
                        (breastCancer_treatment_ind_therapies.study_specific_trt_number=X.study_specific_trt_num5) OR
                        (breastCancer_treatment_ind_therapies.study_specific_trt_number=X.study_specific_trt_num6) ) ) ) as Z
                        
                        ON
                        
                        breastCancer_treatment_dictionary.therapy_chemical_name = Z.therapy_chemical_name ) as Y
                        
                        ON
                        
                        ( (patient_info.study_ID = Y.study_ID) AND  
                        (patient_info.treatment_protocol_number = Y.study_specific_protocol_number) ) 
                        WHERE dbUniquePatientID IN ", dbUniquePatientIDsSQL," ORDER BY patient_info.dbUniquePatientID",sep="")

con<- dbConnect(MySQL(), user="root", dbname="test")
res <- dbSendQuery(con, treatmentQuery)
patientTreatments <- fetch(res, n=-1)
dbDisconnect(con)


patientIDs <- unique(patientTreatments$dbUniquePatientID)


#ASSUMPTION:    the treatmentTable has been populated with your "no" value - for example, 0, NA, etc. 
#I only change "yes" fields.
createPatientTableRow <- function(treatmentTable,GSMID,patientTreatments,classes,drugs,yesField){
  
  dataRows <- patientTreatments[which(patientTreatments$dbUniquePatientID==GSMID),]
  
  tableIndex <- which(rownames(treatmentTable)==GSMID)
  
  
  if( (length(tableIndex)==0) || (length(tableIndex)>1) ){
    
    stop("error! your final treatment table doesn't include this patient ID, or it's duplicated.")
    
  }
  
  for (r in 1:(dim(dataRows)[1])){
    
    matchingTherapyNameIndex <- which(colnames(treatmentTable)==dataRows[r,"therapy_chemical_name"])
    
    if(length(matchingTherapyNameIndex)==0 || length(matchingTherapyNameIndex)>1){
      
      stop("the therapy name in your table pulled from SQL doesn't match up with your treatment table column names or is not specific enough and returned multiple names.")
      
    }
    
    
    #replace the therapy index with a yes field if we haven't encountered this treatment yet for this patient - mark down that yes, this patient took this drug.
    treatmentTable[tableIndex,matchingTherapyNameIndex] <- yesField     
    
    #now look at classes-not just one column per row r!
    #so what are the classes in this data row?
    temp <- match(classes,colnames(dataRows))
    classIndices <- temp[which(!is.na(temp))]
    dataRowsClassesOnly <- dataRows[r,classIndices]
    
    #what classes are actually yes for this patient? "no" fields have already been filled in when we created the treatment table and passed it in.
    yesClassIndices <- which(dataRowsClassesOnly==1)
    #and what are the corresponding indices in the treatment table we're filling in?
    cols <- colnames(dataRowsClassesOnly[r,yesClassIndices])
    temp <- match(cols,colnames(treatmentTable))
    yesclassIndices <- temp[which(!is.na(temp))]
    
    treatmentTable[tableIndex,yesClassIndices] <- yesField
    
    #now add in adjuvant:  if consistent, keep the same. if not, say "mixed"
    #for now, I let an NA be a soft empty case - so if one of the drugs wasn't specified, skip that one.
    #final label with be a drug that WAS actually labeled neo or adj, or mixed.
    if( (treatmentTable[tableIndex,"neoadjuvant_or_adjuvant"]==0) || is.na(treatmentTable[tableIndex,"neoadjuvant_or_adjuvant"]) ){
      #no record yet - probably the first drug we're putting in.
      
      treatmentTable[tableIndex,"neoadjuvant_or_adjuvant"] <- dataRows[r,"neoadjuvant_or_adjuvant"]
      
    }else if(treatmentTable[tableIndex,"neoadjuvant_or_adjuvant"] != dataRows[r,"neoadjuvant_or_adjuvant"]){
      #whoops one drug was neo, one adj: so label this patient as mixed.
      #if they're the same: don't change the label.
      
      treatmentTable[tableIndex,"neoadjuvant_or_adjuvant"] <- "mixed"
      
    }else{
      #just leave it alone - neo or adj tag is consistent across different drugs for this patient.
    }
    
    
  }
  
  
  
  #add in the few last fields:
  
  treatmentTable[tableIndex,"study_ID"] <- dataRows[1,"study_ID"]
  treatmentTable[tableIndex,"study_specific_protocol_number"] <- dataRows[1,"study_specific_protocol_number"]
  
  return(treatmentTable)
  
}


#set up your table for the function above.
colNames <- append(classes,drugs)

#I also want to save a neoadjuvant or adjuvant columns, etc. to fill later for each patient.
#remember this is aggregated across all treatments a patient had: so I can't include dose, etc.
colNames <- append(colNames,c("neoadjuvant_or_adjuvant","study_ID","study_specific_protocol_number") )

noField <- 0;
treatmentTable <- matrix(data=noField, ncol = length(colNames),nrow=length(patientIDs),dimnames=list(patientIDs,colNames))


#now run function to populate table.

for(p in 1:length(patientIDs)){
  
  treatmentTable <- createPatientTableRow(treatmentTable,patientIDs[p],patientTreatments,classes,drugs,yesField=1)
  
}

#rownames are unique patient ID.
output <- list(treatmentTable=treatmentTable,drugs=drugs,classes=classes);

output <- list(treatmentTable=treatmentTable,drugs=drugs,classes=classes);

treatmentTable <- cbind(rownames(treatmentTable),treatmentTable);
treatmentTable <- as.data.frame(treatmentTable);
colnames(treatmentTable)[1] <- "dbUniquePatientID";

fullClinicalTable <- merge(masterPhenoData, treatmentTable,by.x="dbUniquePatientID",by.y="dbUniquePatientID",all.x=TRUE,all.y=TRUE);
#remove the duplicated columns
fullClinicalTable  <- fullClinicalTable[,-c(156,155,154)];
#OK good: study numbers match up.
any(as.character(fullClinicalTable[,152])!=as.character(fullClinicalTable[,2]));
#then remove this column
fullClinicalTable  <- fullClinicalTable[,-c(152)];
#remove any 99832,99831 patients - these are METABRIC.
fullClinicalTable <- fullClinicalTable[-which(fullClinicalTable$study_ID.x==99832), ];
fullClinicalTable <- fullClinicalTable[-which(fullClinicalTable$study_ID.x==99831), ];

clinicalTable <- fullClinicalTable

#remove .x
colnames(clinicalTable)[which(colnames(clinicalTable)=="study_ID.x")] <- "study_ID";
#remove some other extraneous rows
clinicalTable <- clinicalTable[ ,-which(colnames(clinicalTable)=="coordinating_GSE_series_GSMID")];
clinicalTable <- clinicalTable[ ,-which(colnames(clinicalTable)=="microarray_outlier")];
clinicalTable <- clinicalTable[ ,-which(colnames(clinicalTable)=="molecular_subtype_outdatedUsePam50")];
clinicalTable <- clinicalTable[ ,-which(colnames(clinicalTable)=="SET_class")];
clinicalTable <- clinicalTable[ ,-which(colnames(clinicalTable)=="ggi_class")];
clinicalTable <- clinicalTable[ ,-which(colnames(clinicalTable)=="ggi_preTrt")];
clinicalTable <- clinicalTable[ ,-which(colnames(clinicalTable)=="metastasis_stage_endPoint")];
clinicalTable <- clinicalTable[ ,-which(colnames(clinicalTable)=="cellularity")];
clinicalTable <- clinicalTable[ ,-which(colnames(clinicalTable)=="ki67_score")];
#radiotherapy and radiotherapyClass are duplicated.
clinicalTable <- clinicalTable[ ,-which(colnames(clinicalTable)=="radiotherapy")];
#no BRCA info for now.
clinicalTable <- clinicalTable[ ,-which(colnames(clinicalTable)=="BRCA2_mutation")];
clinicalTable <- clinicalTable[ ,-which(colnames(clinicalTable)=="BRCA1_mutation")];
clinicalTable <- clinicalTable[ ,-which(colnames(clinicalTable)=="overall_BRCA_mutation")];
colnames(clinicalTable)[which(colnames(clinicalTable)=="original_study_ID")] <- "original_study_patient_ID"


#add in clinical variable definitions
clinVarDef <- read.csv("~/Dropbox/github/curatedBreastCancer/curatedBreastCancer/clinicalTableNames.csv");
rownames(clinVarDef) <- clinVarDef[ ,1];
#tack on definitions 
clinicalData <- list(clinicalTable=clinicalTable,clinicalVarDef=clinVarDef)
save(clinicalData,file="~/Dropbox/github/curatedBreastCancer/curatedBreastCancer/package/curatedBreastData/data/clinicalData.rda",compress="xz")

#create esets from the RData object I downloaded from MySQL server. each index is a study.
load("~/Dropbox/github/curatedBreastCancer/curatedBreastCancer//allStudies_03_11_2015.RData.gzip")
dataList <- data$allDataSets


#uggh...this is a list! can't really undo unless go column by column.
masterPhenoData <- clinicalTable;

expressionSetList <- list()
exprMatrixList <- list()

#35,36 are METABRIC datasets: private data sharing for this one, can't publically release
#(didn't use in published analyses either - the normalization on this datasets tends to make it a larger outlier.)

for(d in 1:34){
  
  exprMatrixList[[d]] <- list();
  exprMatrixList[[d]]$expr <- dataList[[d]]$expr;
  exprMatrixList[[d]]$featureData <- data.frame(dataList[[d]]$probes,dataList[[d]]$keys);
  colnames(exprMatrixList[[d]]$featureData) <- c("probe","gene_symbol")
  
  #COME BACK: add GEO data - put into protocols??
  #need to clean up IDs so matches with pheno data matrix.
  colnames(exprMatrixList[[d]]$expr) <- substring(colnames(exprMatrixList[[d]]$expr),4,15);
  names(exprMatrixList)[d] <- names(dataList)[d]
  
}

curatedBreastDataExprSetList <- createExpressionSetList(exprMatrixList,masterPhenoData=masterPhenoData,patientKey="GEO_GSMID");

#test: dim(curatedBreastDataExprSetList[[1]]@assayData$exprs)

save(curatedBreastDataExprSetList,file="~/Dropbox/github/curatedBreastCancer/curatedBreastCancer/package/curatedBreastData/data/curatedBreastDataExprSetList.rda",compress="xz");


