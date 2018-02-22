library(REIDS)

#Data
load("HTAData__REIDS_Output.RData")
Probes_and_Junctions=read.table(file="hta1.r1.ass",sep="\t",header=TRUE)

#Groups
groupLiver=c(1,2,3)
groupMuscle=c(4,5,6)

ListofGroups=list(groupLiver,groupMuscle)
names(ListofGroups)=c('groupLiver','groupMuscle')

groups=list(group1=groupLiver,group2=groupMuscle)

HTAData_RASA_PSR_LiverVSMuscle_AllPSR=ASExons(Data=HTAData_OutputREIDSModel,Exonthreshold=0,groups=groups,paired=FALSE,significancelevel=NULL,JunctionInfoFile=Probes_and_Junctions)

save(HTAData_RASA_PSR_LiverVSMuscle_AllPSR,file="HTAData_RASA_PSR_LiverVSMuscle_AllPSR.RData")