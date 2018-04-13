# Processing the HTA CEL Files of the RASA paper
# Create folders
# rawData -> HTA_RASA -> hjay: put in CEL Files
# annotationData -> chipTypes -> hjay : put in CDF File

library(REIDS)

DataProcessing(chipType="HTA-2.0",tags="*r",Name="HTAData",ExonSummarization=TRUE,GeneSummarization=TRUE,FIRMA=TRUE,location="HTAData",verbose=TRUE)

load("HTAData.RData")

PivotTransformData(Data=HTAData,GeneID=HTAData$GeneID,ExonID=HTAData$ExonID,location="HTAData/HTAData_Pivot.csv")