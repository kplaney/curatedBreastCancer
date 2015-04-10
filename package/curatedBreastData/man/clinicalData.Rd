\name{clinicalData}
\alias{clinicalData}
\docType{data}
\title{
Clinical Data Table & Variable Definitions
}
\description{
Clinical data for all samples across all studies, and corresponding variable definitions. Rownames are the GEO_GSMID feature, which corresponds to the sample names in the expression object for a certain study. Includes treatment information.
}
\usage{data("clinicalData")}
\format{
  A list with the following two items:
  -clinicalTable:A data frame. Rownames are the GEO_GSMID feature, which corresponds to the sample names in the expression object for a certain study.  
 -clinicalVarDef:Character string descriptions of each variable.
}
\value{
No return value as this is not a function but rather a data object.
}
\details{
GEO study ID can be found form the study_ID variable. If site_ID is NA, it pertains to the batch ID, which may be due to different platforms being used in the same study or different tissue site collections. Columns 112-151 pertaint to treatment information. radiotherapyClass,  chemotherapyClass, and hormone_therapyClass are indicator variables used to signal whether a patient had radiotherapy, chemotherapy, and/or some form of hormone therapy (usually an estrogen or aromatse inhibitor.) 

More granular information, when available, is provided: for example, whether the chemotherapy drug was capecitabine is coded as the indicator "capecitabine" variable. A value of 1 = yes, 0 = no, NA = not recorded/could not infer from publically available information. "Other" means that most likely, gleaned from the study's Pubmed publication, that the patient may have had other treatments that were not recorded (oftentimes radiotherapy, as this is not always recorded and up to a clinician's discretion in a clinical trial.)

Survival information, such as DFS, RFS, OS, and treatment response information, such pCR and RCB, is also recorded when available.
}

\references{
Planey, Butte. Database integration of 4923 publicly-available samples of breast cancer molecular and clinical data. AMIA Joint Summits Translational Science Proceedings. (2003) PMC3814460
}
\examples{
data(clinicalData)
#check out some of the variable name/definitions
clinicalData$clinicalVarDef[c(1:2),]
#Check out the treatment information. 
#look at first three patients
head(clinicalData$clinicalTable)[c(1:3),c(112:ncol(clinicalData$clinicalTable))]
#how many had chemotherapy?
numChemoPatients <- length(which(
clinicalData$clinicalTable$chemotherapyClass==1))
#how many patients have non-NA OS binary data?
length(which(!is.na(clinicalData$clinicalTable$OS)))
#how many have OS data in the more granular form of months until OS? 
#this variable includes studies that had a cieling for tracking OS
length(which(!is.na(clinicalData$clinicalTable$OS_months_or_MIN_months_of_OS)))
#how many patients have OS information that is definitively 
#followed up until their death
#(details on how studies collect OS data can be surprising!)
length(which(!is.na(clinicalData$clinicalTable$OS_up_until_death)))

}
\keyword{datasets}
