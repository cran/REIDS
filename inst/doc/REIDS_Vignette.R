### R code from vignette source 'REIDS_Vignette.Rnw'

###################################################
### code chunk number 1: DataProcessing (eval = FALSE)
###################################################
## DataProcessing(chipType="HuEx-1_0-st-v2",tags="coreR3,A20071112,EP",
##                Name="TissueData",ExonSummarization=TRUE,GeneSummarization=TRUE,
##                FIRMA=TRUE,location="TissueData",verbose=TRUE)


###################################################
### code chunk number 2: PivotTransformation (eval = FALSE)
###################################################
## load("TissueData/TissueData.RData")
## 
## PivotTransformation(Data=TissueData,GeneID=TissueData$GeneID,
##                     ExonID=TissueData$ExonID,savecsv=TRUE,
##                     Name="TissueData_Pivot",location="TissueData")


###################################################
### code chunk number 3: REIDS (eval = FALSE)
###################################################
## load("TissueData/TissueData.RData")
## 
## REIDS_Output=REIDSFunction(geneData=TissueData,nsim=5000,geneID=TissueData$GeneID,
##                        exonID=TissueData$ExonID,informativeCalls=TRUE,alpha=0.5)


###################################################
### code chunk number 4: Search (eval = FALSE)
###################################################
## exonID <- c(3252129,3597384,3333718,3735208,2598321,3338589,2605391,2605390,
##             2605386,3025632,2375766,3569827,3569830,2334499,3972987,2516011,
##             2989068,3422189)
## 
## load("REIDS_Output.RData")
## 
## Results=Search(WhatToLookFor=data.frame(ExonID=exonID),Data=REIDS_Output, 
##                AggregateResults=FALSE,NotFound=NULL)


###################################################
### code chunk number 5: ExonTesting (eval = FALSE)
###################################################
## load("TissueData/REIDS_Output.RData")
## 
## groupHMPT=c(7,8,9,16,17,18,22,23,24,31,32,33)
## groupOthers=c(1,2,3,4,5,6,10,11,12,13,14,15,19,20,21,25,26,27,28,29,30)
## groups=list(group1=groupHMPT,group2=groupOthers)
## 
## TissueData_ExonTesting=ExonTesting(Data=REIDS_Output,Exonthreshold=NULL,groups=groups,
##                                    paired=FALSE,significancelevel=NULL)


###################################################
### code chunk number 6: ASExons (eval = FALSE)
###################################################
## TissueData_ExonTesting_Sign1=ASExons(Data=REIDS_Output,Exonthreshold=0.5,groups=groups,
##                                      paired=FALSE,significancelevel=0.05,JunctionFile=NULL)
## 
## TissueData_ExonTesting_Sign2=ASExons(Data=TissueData_ExonTesting,Exonthreshold=0.5,
##                                      groups=groups,paired=FALSE,significancelevel=0.05
## 									 ,JunctionFile=NULL)


###################################################
### code chunk number 7: FIRMAScores (eval = FALSE)
###################################################
## load("TissueData/FIRMA_Output.RData")
## load("TissueData/REIDS_Output.RData")
## 
## groupHMPT=c(7,8,9,16,17,18,22,23,24,31,32,33)
## groupOthers=c(1,2,3,4,5,6,10,11,12,13,14,15,19,20,21,25,26,27,28,29,30)
## groups=list(group1=groupHMPT,group2=groupOthers)
## 
## exons=TissueData_ExonTesting[which(TissueData_ExonTesting$X50.>0.5),]
## 
## TissueData_FIRMA=FIRMAScores(Data=FIRMA_Output,InformativeExons=exons,
##                              groups=groups,paired=FALSE,significancelevel=0.05)


###################################################
### code chunk number 8: PlotFunction (eval = FALSE)
###################################################
## #Top gene 2736322 and   Top exon 2736397
## load("TissueData/TissueData.rda")
## load("TissueData/TissueData_ExonLevelSummarized.RData")
## load("TissueData/TissueData_GeneLevelSummarized.RData")
## load("TissueData/FIRMA_Output.RData")
## load("TissueData/REIDS_Output.RData") 
## 
## groupHMPT=c(7,8,9,16,17,18,22,23,24,31,32,33)
## groupOthers=c(1,2,3,4,5,6,10,11,12,13,14,15,19,20,21,25,26,27,28,29,30)
## groups=list(group1=groupHMPT,group2=groupOthers)
## 
## PlotFunction(GeneID="2736322",ExonID="2736397",Data=TissueData,
##              REIDS_Output=REIDS_Output,GeneLevelData=TissueData_GeneLevelSummarized,
##              ExonLevelData=TissueData_ExonLevelSummarized,FIRMA_Data=FIRMA_Output,
##              groups=groups,plottype="new",
##              location="TissueData/Gene2736322_Exon2736397")
## 


###################################################
### code chunk number 9: SpliceIndex (eval = FALSE)
###################################################
## load("TissueData/TissueData_ExonLevelSummarized.RData")
## load("TissueData/TissueData_GeneLevelSummarized.RData")
## 
## 
## groupHMPT=c(7,8,9,16,17,18,22,23,24,31,32,33)
## groupOthers=c(1,2,3,4,5,6,10,11,12,13,14,15,19,20,21,25,26,27,28,29,30)
## groups=list(group1=groupHMPT,group2=groupOthers)
## 
## 
## SI_Output=SpliceIndex(GeneData=TissueData_ExonLevelSummarized,
##                       ExonData=TissueData_GeneLevelSummarized,
##                       InformativeExons=NULL,groups=groups,
##                       paired=FALSE,significancelevel=NULL)


###################################################
### code chunk number 10: HTAData (eval = FALSE)
###################################################
## library(REIDS)
## 
## DataProcessing(chipType="hjay",tags="coreR3,A20071112,EP",Name="HTAData",
## 		ExonSummarization=TRUE,GeneSummarization=TRUE,FIRMA=TRUE,
## 		location="HTAData",verbose=TRUE)
## 
## load("HTAData.RData")
## 
## PivotTransformData(Data=HTAData,GeneID=HTAData$GeneID,ExonID=HTAData$ExonID,
## 		savecsv=FALSE,Name="HTAData_Pivot",location="HTAData")


###################################################
### code chunk number 11: REIDSHTA (eval = FALSE)
###################################################
## load("TissueData/TissueData.RData")
## 
## REIDS_Output=REIDSFunction(geneData=TissueData,nsim=5000,geneID=TissueData$GeneID,
##                        exonID=TissueData$ExonID,informativeCalls=TRUE,alpha=0.5)


###################################################
### code chunk number 12: HTA_Analysis (eval = FALSE)
###################################################
## load("HTAData__REIDS_Output.RData")
## Probes_and_Junctions=read.table(file="hta1.r1.ass",sep="\t",header=TRUE)
## 
## groupLiver=c(1,2,3)
## groupMuscle=c(4,5,6)
## 
## ListofGroups=list(groupLiver,groupMuscle)
## names(ListofGroups)=c('groupLiver','groupMuscle')
## 
## HTAData_LiverVSMuscle=ASExons(Data=HTAData_OutputREIDSModel,Exonthreshold=0,
## 		groups=groups,paired=FALSE,significancelevel=NULL,
## 		JunctionInfoFile=Probes_and_Junctions)
## 
## save(HTAData_LiverVSMuscle,file="HTAData_LiverVSMuscle.RData")


###################################################
### code chunk number 13: HTARank (eval = FALSE)
###################################################
## load("Data/HTAData_LiverVSMuscle.RData")
## load("Data/HTAData_RASA.rda")
## for(i in 3:ncol(HTAData_RASA)){
## 	HTAData_RASA[,i]=as.numeric(as.character(HTAData_RASA[,i]))
## }
## 
## Filter=HTAData_LiverVSMuscle[HTAData_LiverVSMuscle$X50.>=0.50,]
## 
## ASPSR_PSR=Filter[Filter$adj.p.value<0.05,]
## length(unique(ASPSR_PSR$ExonID))
## 
## REIDS_RankedProbesets=REIDS_JunctionAssessment(DABGFile,probeset_probesfile,
## 					ASProbeSets,AnnotData,Data,mode=c("Liberal"))
## save(REIDS_RankedProbesets,file="REIDS_RankedProbesets.RData")


###################################################
### code chunk number 14: sessionInfo
###################################################
toLatex(sessionInfo())


