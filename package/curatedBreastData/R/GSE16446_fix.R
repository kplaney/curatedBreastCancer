#the fix I made when a user noted that for GSE16446, the DFS and OS values were still in days, not months.
#you can confirm this on the GEO website - the original values were in days.
#this change was pushed to Bioconductor on 03/04/17 and I will watch to make sure it's incorporated in the 
#next release

#assumes you are in curatedBreastData/data/

load("curatedBreastDataExprSetList.rda")
#fixing GSE16446
tmp = pData(curatedBreastDataExprSetList[["study_16446_GPL570_all"]])
#need to go from days to months, so divide by 28
tmp$DFS_months_or_MIN_months_of_DFS = tmp$DFS_months_or_MIN_months_of_DFS/28
tmp$OS_months_or_MIN_months_of_OS = tmp$OS_months_or_MIN_months_of_OS/28
pData(curatedBreastDataExprSetList[["study_16446_GPL570_all"]]) = tmp
validObject(curatedBreastDataExprSetList[["study_16446_GPL570_all"]], "eset")
save(curatedBreastDataExprSetList, file = "curatedBreastDataExprSetList.rda", compress = "xz")
rm(curatedBreastDataExprSetList)


load("clinicalData.rda")
clinicalData$clinicalTable[clinicalData$clinicalTable$study_ID == 16446, ]$DFS_months_or_MIN_months_of_DFS = 
  clinicalData$clinicalTable[clinicalData$clinicalTable$study_ID == 16446, ]$DFS_months_or_MIN_months_of_DFS/28

clinicalData$clinicalTable[clinicalData$clinicalTable$study_ID == 16446, ]$OS_months_or_MIN_months_of_OS = 
  clinicalData$clinicalTable[clinicalData$clinicalTable$study_ID == 16446, ]$OS_months_or_MIN_months_of_OS/28
save(clinicalData, file = "clinicalData.rda", compress = "xz")
rm(clinicalData)
#make sure to bump the version of the file, too
#changed to Version: 2.3.0 in DESCRIPTION file.
