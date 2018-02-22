# Processing the HTA CEL Files of the RASA paper
# Create folders
# rawData -> HTA_RASA -> hjay: put in CEL Files
# annotationData -> chipTypes -> hjay : put in CDF File

library(REIDS)

DataProcessing(chipType="hjay",tags="coreR3,A20071112,EP",Name="HTAData",ExonSummarization=TRUE,GeneSummarization=TRUE,FIRMA=TRUE,location="HTAData",verbose=TRUE)

load("HTAData.RData")

PivotTransformData(Data=HTAData,GeneID=HTAData$GeneID,ExonID=HTAData$ExonID,savecsv=FALSE,Name="HTAData_Pivot",location="HTAData")