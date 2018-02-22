library(REIDS)

load("HTAData_ExonLevelSummarized.RData")
load("HTAData_GeneLevelSummarized.RData")


groupLiver=c(1,2,3)
groupMuscle=c(4,5,6)
ListofGroups=list(groupLiver,groupMuscle)
names(ListofGroups)=c('groupLiver','groupMuscle')

groups=list(group1=groupLiver,group2=groupMuscle)
HTAData_SI_Output=SpliceIndex(GeneData=HTAData_ExonLevelSummarized,ExonData=HTAData_GeneLevelSummarized,InformativeExons=NULL,groups=groups,paired=FALSE,significancelevel=NULL)
save(HTAData_SI_Output,file="HTAData_SI_Output.RData")