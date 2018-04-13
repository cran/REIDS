library(REIDS)

#Data
load("HTAData__REIDS_Output.RData")
Probes_and_Junctions=read.table(file="hta1.r1.ass",sep="\t",header=TRUE)

#Groups
Groups=list(c(10:18),c(19:27))

HTAData_1vs2=ASExons(Data=HTA-2.0,REIDS_Output,Exonthreshold=0.5,Groups=list(c(10:18),c(19:27)),paired=FALSE,significancelevel=0.05)

save(HTAData_1vs2,file="HTAData_1vs2.RData")