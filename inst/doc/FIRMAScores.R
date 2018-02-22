library(REIDS)

load("FIRMA_Output")

groupLiver=c(1,2,3)
groupMuscle=c(4,5,6)
ListofGroups=list(groupLiver,groupMuscle)
names(ListofGroups)=c('groupLiver','groupMuscle')
groups=list(group1=groupLiver,group2=groupMuscle)


FIRMA_Output_Scores=FIRMAScores(Data=FIRMA_Output,InformativeExons=NULL,groups=groups,paired=FALSE,significancelevel=NULL)

save(FIRMA_Output_Scores,file="FIRMA_Output_Scores.RData")