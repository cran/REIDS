### R code from vignette source 'REIDS_Vignette.Rnw'

###################################################
### code chunk number 1: DataProcessing (eval = FALSE)
###################################################
## library(REIDS)
## DataProcessing(chipType="HTA-2_0",tags="*,r",Name="HTAData",
## 		ExonSummarization=TRUE,GeneSummarization=TRUE,FIRMA=TRUE,
## 		location="HTAData",verbose=TRUE)


###################################################
### code chunk number 2: PivotTransformation (eval = FALSE)
###################################################
## load("HTAData.RData")
## 
## PivotTransformData(Data=HTAData,GeneID=HTAData$GeneID,ExonID=HTAData$ExonID,
## 		savecsv=FALSE,Name="HTAData_Pivot",location="HTAData")
## 
## GeneID=unique(as.character(HTAData$GeneID))
## GeneTable=data.frame("GeneID"=GeneID)
## save(GeneTable,file="GeneTable.RData")


###################################################
### code chunk number 3: DABGProbesets (eval = FALSE)
###################################################
## ProbesetID_ProbesetName=read.table("PSID_Name_TCID.txt")
## DABG=read.table("APT_DABG/dabg.summary.txt",header=TRUE,sep="\t")
## 
## #DABG in all samples
## Low=apply(DABG,1,function(x) all(x>=0.05))
## LowUIDs=DABG[,1][Low]
## LOW_Probesets_All=unique(ProbesetID_ProbesetName[which(
## ProbesetID_ProbesetName[,1]%in%LOWUIDs),2])
## save(LOW_Probesets_All,file="Low_AllSamples.RData")
## 
## #Per Group
## Low_GSamples=list()
## for(g in 1:length(Groups)){
## 	Low=apply(DABG,1,function(x) all(x[Groups[[g]]]>=0.05))
## 	LowUIDs=DABG[,1][Low]
## 	LOW_Probesets_Samples=unique(ProbesetID_ProbesetName[which(
## 	ProbesetID_ProbesetName[,1]%in%LOWUIDs),2])
## 	Low_GSamples[[g]]=LOW_Probesets_Samples
## }
## save(Low_GSamples,file="Low_GSamples.RData")


###################################################
### code chunk number 4: REIDS (eval = FALSE)
###################################################
## load("GeneTable.RData")
## load("Low_AllSamples.RData")
## 
## 
## REIDSFunction(geneIDs=GeneTable,Name="HTA-2.0", Indices="HTAData_lineindex.csv",
## 		DataFile="HTAData_Pivot.csv",nsim=5000,informativeCalls=FALSE,
## 		rho=0.5,Low_AllSamples=Low_AllSamples,Location="Output")
## 
## 


###################################################
### code chunk number 5: Search (eval = FALSE)
###################################################
## exonID <- c("PSR01003414","PSR02012920","PSR12000150")
## 
## load("HTA-2.0_REIDS_Output.RData")
## 
## Results=Search(WhatToLookFor=data.frame(ExonID=exonID),Data=REIDS_Output, 
##                AggregateResults=FALSE,NotFound=NULL)


###################################################
### code chunk number 6: ExonTesting (eval = FALSE)
###################################################
## load("HTA-2.0_REIDS_Output.RData")
## 
## HTAData_1vs2=ASExons(Data=HTA-2.0,REIDS_Output,Exonthreshold=0.5,
## 		Groups=list(c(10:18),c(19:27)),paired=FALSE,
## 		significancelevel=0.05)
## 
## save(HTAData_1vs2,file="HTAData_1vs2.RData")


###################################################
### code chunk number 7: HTAJUN (eval = FALSE)
###################################################
## load("HTAData_1vs2.RData")
## load("HTAData.RData")
## load("Low_GSamples.RData")
## load("Low_AllSamples.RData")
## for(i in 3:ncol(HTAData)){
## 	HTAData[,i]=as.numeric(as.character(HTAData[,i]))
## }
## 
## ASPSR_PSR=HTAData_1vs2[HTAData_1vs2$adj.p.value<0.05,]
## length(unique(ASPSR_PSR$ExonID))
## 
## #ExonAnnotations
## EI=read.table("LineIndexing_ExonAnnot.txt",header=FALSE,sep="\t",
## 		stringsAsFactors=FALSE)
## 
## #JunctionAnnotations
## JI=read.table("LineIndexing_JAnnot.txt",header=FALSE,sep="\t",
## 		stringsAsFactors=FALSE)
## 
## #Transcript Annotations
## TrI=read.table("LineIndexing_TrAnnot.txt",header=FALSE,sep="\t",
## 		stringsAsFactors=FALSE)
## 
## 
## REIDS_JunctionAssesment(ASProbeSets=unique(ASPSR_PSR$ExonID),
## 		JAnnotI=JI,EAnnotI=EI,TrAnnotI=TrI,Data=HTAData,Groups=list(c(3,4,5),
## 		c(6,7,8)),Low_AllSamples,Low_GSamples,Plot=FALSE,Name="HTA-2.0")
## 	


###################################################
### code chunk number 8: PlotFunction1 (eval = FALSE)
###################################################
## load("HTAData.RData")
## load("ExonLevel.RData")
## load("GeneLevel.RData")
## ExpressionLevelPlot(GeneID="TC12000010",ExonID="PSR12000150",Data=HTAData,
## 		GeneLevelData=GeneLevel,ExonLevelData=ExonLevel,groups=
## 		list(c(1:3),c(4:6)),ylabel="",title="PSR12000150")	
## 


###################################################
### code chunk number 9: PlotFunction2 (eval = FALSE)
###################################################
## load("HTAData.RData")
## load("ExonLevel.RData")
## 
## library(data.table) 
## transcript.clusters.NetAffx.36<-fread("HTA-2_0.na35.2.hg19.transcript.csv",
## 		skip=19,header=TRUE,sep=",")
## transcript.clusters.NetAffx.36=as.data.frame(transcript.clusters.NetAffx.36)
## transcript.clusters.NetAffx.36=transcript.clusters.NetAffx.36[,c(1,4,8)]
## save(transcript.clusters.NetAffx.36,file="transcript.clusters.NetAffx.36.RData")
## #load("transcript.clusters.NetAffx.36.RData")
## 
## probesets.NetAffx.36 <- fread("HTA-2_0.na35.hg19.probeset.csv",
## 		skip=19,header=TRUE,sep=",") 
## probesets.NetAffx.36 =as.data.frame(probesets.NetAffx.36)
## positions_36 <- probesets.NetAffx.36[probesets.NetAffx.36[,13] == "main",c(1,2,4,5,6,7)]
## save(positions_36,file="positions_36.RData")
## #load("positions_36.RData")
## 
## 
## TranscriptsPlot(trans="TC12000010", positions=positions_36,
## 		transcriptinfo=transcript.clusters.NetAffx.36,
## 		display.probesets=TRUE,Data=ExonLevel,groups=list(c(1:3),c(4:6)),
## 		Start=NULL,Stop=NULL,Highlight="PSR12000150")	
## 


###################################################
### code chunk number 10: PlotFunction3 (eval = FALSE)
###################################################
## load("Data/HTAData_1vs2.RData")
## load("Data/HTAData.RData")
## load("Low_GSamples.RData")
## load("Low_AllSamples.RData")
## for(i in 3:ncol(HTAData_RASA)){
## 	HTAData_RASA[,i]=as.numeric(as.character(HTAData_RASA[,i]))
## }
## 
## ExonScoreFilter=HTAData_LiverVSMuscle[HTAData_LiverVSMuscle$X50.>=0.50,]
## 
## ASPSR_PSR=ExonScoreFilter[ExonScoreFilter$adj.p.value<0.05,]
## length(unique(ASPSR_PSR$ExonID))
## 
## #ExonAnnotations
## EI=read.table("LineIndexing_ExonAnnot.txt",header=FALSE,sep="\t",
## 		stringsAsFactors=FALSE)
## 
## #JunctionAnnotations
## JI=read.table("LineIndexing_JAnnot.txt",header=FALSE,sep="\t",
## 		stringsAsFactors=FALSE)
## 
## #Transcript Annotations
## TrI=read.table("LineIndexing_TrAnnot.txt",header=FALSE,sep="\t",
## 		stringsAsFactors=FALSE)
## 
## JunInfo(x="TC12000010",ASPSR=unique(ASPSR_PSR$ExonID),JLines=JI[which(JI[,1]==x),],
## 		TrLines=TrI[which(TrI[,1]==x),],ELines=I[which(EI[,1]==x),],
## 		DataS=HTADataData[which(as.character(HTAData[,1])==x),],
## 		Groups=list(c(1:3),c(4:6)),Low_ALLSamples=Low_ALLSamples,Low_GSamples=Low_GSamples,Plot=TRUE)


###################################################
### code chunk number 11: sessionInfo
###################################################
toLatex(sessionInfo())


